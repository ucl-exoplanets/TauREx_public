################################################
#class fitting
#
# Class containing fitting interfaces to pymc and multinest.
# All likelihoods are defined here.  
#
# Input: -parameter object
#        -data object
#        -
#
#
# Output: -data object containing relevant data in dictionary 
#         format  
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013 
#
################################################

#loading libraries     
from base import base
import numpy, pylab, os, sys, math, pymc, warnings, threading, subprocess, gzip, pickle, shutil, logging
from pylab import *
from numpy import *
from StringIO import StringIO
from scipy.interpolate import interp1d
from scipy.optimize import fmin
from scipy.optimize import minimize

try: 
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

try:
    from mpi4py import MPI
    MPIimport = True
except ImportError:
    MPIimport = False
    pass

class fitting(base):

    def __init__(self, forwardmodel, params=None, data=None, atmosphere=None, rad_model=None):

        # note that forwardmodel can be either the transmission or the emission object

        self.__ID__ = 'fitting'

        logging.info('Initialise object fitting')

        if params:
            self.params = params
        else:
            self.params = forwardmodel.params # get params object from profile

        if data:
            self.data = data
        else:
            self.data = forwardmodel.data    # get data object from forwardmodel

        if atmosphere:
            self.atmosphere = atmosphere
        else:
            self.atmosphere = forwardmodel.atmosphere    # get atmosphere object from forwardmodel

        self.forwardmodel = forwardmodel

        self.forwardmodel_type = type(forwardmodel).__name__
        logging.info('Radiative transfer model: %s' % self.forwardmodel_type)

        # MPI support
        if MPIimport:
            self.MPIrank     = MPI.COMM_WORLD.Get_rank()
            self.MPIsize     = MPI.COMM_WORLD.Get_size()
        else:
            self.MPIrank     = 0
            self.MPIsize     = 0

        # set some folder names
        self.dir_mcmc = os.path.join(self.params.out_path, 'MCMC')
        self.dir_multinest = os.path.join(self.params.out_path, 'multinest')

        # create some folders.
        if self.MPIrank == 0:
            folders = [self.params.out_path, self.dir_mcmc, self.dir_multinest]
            for f in folders:
                if not os.path.isdir(f):
                    logging.info('Create folder %s' % f)
                    os.mkdir(f)

        #getting parameter grid with initial estimates
        self.fit_params_init = self.atmosphere.fit_params
        self.fit_index = self.atmosphere.fit_index
        self.fit_nparams = len(self.atmosphere.fit_params)

        #getting prior bounds for downhill algorithm
        self.bounds = self.atmosphere.bounds

        #initialising output tags
        self.DOWNHILL = False
        self.MCMC = False
        self.NEST = False

        #DEBUG COUNTER
        self.db_count = 0

    def chisq_trans(self, fit_params, data, datastd):

        #chisquare minimisation bit

        T, P, X = self.atmosphere.TP_profile(fit_params=fit_params)  # calculating TP profile

        # calculating densities
        rho = self.atmosphere.get_rho(T=T, P=P)

        #the temperature parameter should work out of the box but check for transmission again
        model = self.forwardmodel.model(rho=rho, X=X, temperature=T)

        #binning internal model
        model_binned = [model[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid)]
        
        ion()
        figure(1)
        clf()
        plot(self.data.spectrum[:,0],self.data.spectrum[:,1])
        plot(self.data.spectrum[:,0], model_binned)
        draw()

        res = (data - model_binned) / datastd

        return sum(res**2)

###############################################################################
#simplex downhill algorithm

    def downhill_fit(self):

        logging.info('Fit data using %s minimisation' % self.params.downhill_type)

        # @todo here we could use different starting points for different threads?

        data = self.data.spectrum[:,1] # observed data
        datastd = self.data.spectrum[:,2] # data error

        fit_output = minimize(fun=self.chisq_trans,
                              x0=self.fit_params_init,
                              args=(data,datastd),
                              method=self.params.downhill_type,
                              bounds=(self.bounds))

        Tout_mean, Xout_mean = self.collate_downhill_results(fit_output['x'])

        self.DOWNHILL = True
        self.DOWNHILL_T_mean = Tout_mean
        self.DOWNHILL_X_mean = Xout_mean
        self.DOWNHILL_fit_output = fit_output['x']
        
        print self.DOWNHILL_fit_output

    def collate_downhill_results(self, fit_params):

        logging.info('Unpacking the downhill minimization results')

        T, P, X = self.atmosphere.TP_profile(fit_params=fit_params)

        return T, X


###############################################################################    
#Markov Chain Monte Carlo algorithm


    # noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment
    # @profile #line-by-line profiling decorator
    def mcmc_fit(self):
    # adaptive Markov Chain Monte Carlo

        logging.info('Start MCMC fit')

        if self.DOWNHILL:
            fit_params_init = self.DOWNHILL_fit_output
        else:
            fit_params_init = self.fit_params_init

        data = self.data.spectrum[:,1] #observed data
        datastd = self.data.spectrum[:,2] #data error
        
        # setting prior distributions
        # master thread (MPIrank =0) will start from ideal solution (from downhill fitting)
        # slave threads (MPIrank != 0) will start from randomised points.

        priors = empty(size(fit_params_init), dtype=object)

        #setting up main thread. Use downhill FIT as starting points
        if self.MPIrank == 0:
            for i in range(self.fit_nparams):
                    priors[i] = pymc.Uniform('PFIT_%i' % (i), self.bounds[i][0], self.bounds[i][1],
                                             value=fit_params_init[i])  # uniform prior

        else:
            #setting up other threads (if exist). Their initial starting positions will be randomly perturbed
            for i in range(self.fit_nparams):
                param_range = (self.bounds[i][1] - self.bounds[i][0]) / 5.0 #range of parameter over which to perturb starting position
                param_mean  = np.mean(self.bounds[i])
                param_rand  = random.uniform(low=param_mean-param_range,high=param_mean+param_range) #random parameter start
                priors[i] = pymc.Uniform('PFIT_%i' % (i), self.bounds[i][0],self.bounds[i][1],value=param_rand)  # uniform prior

        #setting up data error prior if specified
        if self.params.mcmc_update_std:
            #uniform prior on data standard deviation
            std_dev = pymc.Uniform('datastd', 0.0, 2.0*max(datastd), value=datastd, size=len(datastd))
        else:
            std_dev = pymc.Uniform('datastd', 0.0, 2.0*max(datastd), value=datastd, observed=True, size=len(datastd))


        def mcmc_loglikelihood(value, fit_params, datastd, data):
            # log-likelihood function. Needs to be initialised directly since CYTHON does not like PYMC decorators
            # @todo need to cq ast from numpy object array to float array. Slow but hstack is  slower.

            fit_params_container = np.zeros((self.fit_nparams))
            for i in range(self.fit_nparams):
                fit_params_container[i] = fit_params[i]
            # @todo params_container should be equal to PFIT? I think so... Maybe not if we fix some values

            chi_t = self.chisq_trans(fit_params_container, data, datastd) #calculate chisq
            llterms =   (-len(data)/2.0)*np.log(pi) - np.log(np.mean(datastd)) - 0.5* chi_t

            return llterms

        mcmc_logp = pymc.Stochastic(logp = mcmc_loglikelihood,
                                    doc = 'The switchpoint for mcmc loglikelihood.',
                                    name = 'switchpoint',
                                    parents = {'fit_params': priors, 'datastd': datastd, 'data': data},
                                    #random = switchpoint_rand,
                                    trace = True,
                                    value = fit_params_init,
                                    dtype=int,
                                    rseed = 1.,
                                    observed = True,
                                    cache_depth = 2,
                                    plot=False,
                                    verbose = 0)

        # set output folder
        dir_mcmc_thread= os.path.join(self.dir_mcmc, 'thread_%i' % self.MPIrank)

        # remove old files
        if os.path.isdir(dir_mcmc_thread):
            logging.debug('Remove folder %s' % dir_mcmc_thread)
            shutil.rmtree(dir_mcmc_thread)

        # executing MCMC sampling
        if self.params.mcmc_verbose:
            verbose = 1
        else:
            verbose = 0

        # build the model
        R = pymc.MCMC((priors, mcmc_logp),
                      verbose=verbose,
                      db='txt',
                      dbname=dir_mcmc_thread)

        # populate and run it
        R.sample(iter=self.params.mcmc_iter,
                 burn=self.params.mcmc_burn,
                 thin=self.params.mcmc_thin,
                 progress_bar=self.params.mcmc_progressbar,
                 verbose=verbose)

        #coallating results into arrays
        Tout_mean, Tout_std, Xout_mean, Xout_std = self.collate_mcmc_result(R)

        #saving arrays to object. This should really be divided by rank, see below
        self.MCMC = True
        self.MCMC_T_mean = Tout_mean
        self.MCMC_T_std = Tout_std
        self.MCMC_X_mean = Xout_mean
        self.MCMC_X_std = Xout_std
        self.MCMC_STATS = R.stats()
        self.MCMC_FITDATA= R

        # each MCMC chain is run on a separate thread. Save all outputs
        logging.info('Store the MCMC results')
        self.MCMC_OUT = {}
        self.MCMC_OUT[self.MPIrank] = {}
        self.MCMC_OUT[self.MPIrank]['FITDATA'] = R
        self.MCMC_OUT[self.MPIrank]['STATS'] = R.stats()
        self.MCMC_OUT[self.MPIrank]['T_mean'] = Tout_mean
        self.MCMC_OUT[self.MPIrank]['T_std'] = Tout_std
        self.MCMC_OUT[self.MPIrank]['X_mean'] = Xout_mean
        self.MCMC_OUT[self.MPIrank]['X_std'] = Xout_std

    def collate_mcmc_result(self, MCMCout):

        logging.info('Unpacking the MCMC results')

        MCMCstats = MCMCout.stats()

        fit_params_out = []
        fit_params_out_std = []
        for i in range(len(self.fit_params_init)):
            fit_params_out.append(MCMCstats['PFIT_%i' % i]['mean'])
            fit_params_out_std.append(MCMCstats['PFIT_%i' % i]['standard deviation'])

        T,P,X = self.atmosphere.TP_profile(fit_params=fit_params_out)

        T_std = fit_params_out_std[self.fit_index[0]:self.fit_index[1]]
        X_std = fit_params_out_std[:self.fit_index[0]]

        return np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std)


############################################################################### 
#Nested Sampling algorithm

    def multinest_fit(self, resume=None):
        #multinest fitting routine (wrapper around PyMultiNest)

        logging.info('Start MULTINEST fit')

        if resume is None:
            resume = self.params.nest_resume

        if resume:
            filename = os.path.join(self.dir_multinest, '1-live.points')
            if os.path.isfile(filename):
                livepoints = sum(1 for line in open(filename))
                if livepoints != self.params.nest_nlive:
                    logging.warning('Cannot resume previous MULTINEST run, the number of live points has changed')
                    logging.warning('Removing previous MULTINEST chains')
                    for file in os.listdir(self.dir_multinest):
                        file_path = os.path.join(self.dir_multinest, file)
                        try:
                            if os.path.isfile(file_path):
                                os.unlink(file_path)
                        except Exception, e:
                            logging.error('Cannot remove files in %s' % self.dir_multinest)

        def show(filepath): 
            # open the output (pdf) file for the user
            if os.name == 'mac':
                subprocess.call(('open', filepath))
            elif os.name == 'nt':
                os.startfile(filepath)

        def multinest_loglike(cube, ndim, nparams):
            #log-likelihood function called by multinest
            fit_params_container = asarray([cube[i] for i in xrange(self.fit_nparams)])
            chi_t = self.chisq_trans(fit_params_container, data, datastd)
            llterms = (-ndim/2.0)*np.log(2.*pi*datastd_mean**2) - 0.5*chi_t
            return llterms
    
        def multinest_uniform_prior(cube, ndim, nparams):
            #prior distributions called by multinest. Implements a uniform prior
            #converting parameters from normalised grid to uniform prior
            for i in xrange(self.fit_nparams):
                cube[i] = (cube[i] * (self.bounds[i][1]-self.bounds[i][0])) + self.bounds[i][0]

        data = self.data.spectrum[:,1] #observed data
        datastd = self.data.spectrum[:,2] #data error
        datastd_mean = mean(datastd)
        
        parameters = [str(i) for i in xrange(self.fit_nparams)]
        ndim = self.fit_nparams

        #progress = pymultinest.ProgressPlotter(n_params = self.fit_nparams); progress.start()
        #threading.Timer(60, show, ["chains/1-phys_live.points.pdf"]).start() # delayed opening
        logging.debug('Multinest output dir %s' % self.dir_multinest)

        pymultinest.run(LogLikelihood=multinest_loglike,
                        Prior=multinest_uniform_prior,
                        n_dims=self.fit_nparams,
                        multimodal=self.params.nest_multimodes,
                        max_modes=self.params.nest_max_modes,
                        outputfiles_basename=os.path.join(self.dir_multinest, '1-'),
                        const_efficiency_mode = self.params.nest_const_eff,
                        importance_nested_sampling = self.params.nest_imp_sampling,
                        resume = resume,
                        verbose = self.params.nest_verbose,
                        sampling_efficiency = self.params.nest_samp_eff,
                        n_live_points = self.params.nest_nlive,
                        max_iter= self.params.nest_max_iter,
                        init_MPI=False)
        #progress.stop()

        # wait for all threads to synchronise
        MPI.COMM_WORLD.Barrier()

        if MPI.COMM_WORLD.Get_rank() == 0:

            #coallating results into arrays (only for the main thread)

            logging.info('Store the MULTINEST results')
            NESTout = pymultinest.Analyzer(n_params=self.fit_nparams,
                                       outputfiles_basename=os.path.join(self.dir_multinest, '1-'))

            # this returns a list of solutions: [(T1, T_std1, X1, X_std1), (T2 T_std2, X2, X_std2), ... ]
            multinest_result = self.collate_multinest_result(NESTout)

            #saving arrays to object
            self.NEST = True
            self.NEST_modes = len(multinest_result) # number of modes detected
            self.NEST_T_mean = [solution[0] for solution in multinest_result]
            self.NEST_T_std = [solution[1] for solution in multinest_result]
            self.NEST_X_mean = [solution[2] for solution in multinest_result]
            self.NEST_X_std = [solution[3] for solution in multinest_result]
            self.NEST_stats = NESTout.get_stats()
            self.NEST_FITDATA = NESTout

    def collate_multinest_result(self, NESTout):

        # @todo to think about!
        # Multinest in multimdes can return two or more sets of solutions.
        # if only one solution is found, it returns a list with a tuple with 4 elements: [(T, T_std, X, X_std)]
        # if N (>1) solutions are found, it returns a list with N tuples with 4 elements each
        #    e.g. [(T1, T_std1, X1, X_std1), (T2 T_std2, X2, X_std2), ... ]

        logging.info('Unpacking the MULTINEST results')

        NESTstats = NESTout.get_stats()

        if self.params.nest_multimodes:
            if len(NESTstats['modes']) > 1:
                # more than one mode found
                out_list = []
                for n in range(len(NESTstats['modes'])):
                    fit_params_out = []
                    fit_params_out_std = []
                    for i in range(self.fit_nparams):
                        fit_params_out.append(NESTstats['modes'][n]['maximum a posterior'][i])
                        fit_params_out_std.append(NESTstats['modes'][n]['sigma'][i])

                    T, P, X = self.atmosphere.TP_profile(fit_params=fit_params_out)
                    T_std = fit_params_out_std[self.fit_index[0]:self.fit_index[1]]
                    X_std = fit_params_out_std[:self.fit_index[0]]
                    out_list.append((np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std)))

                return out_list

        # only one set of solutions found
        fit_params_out = []
        fit_params_out_std = []
        for i in range(self.fit_nparams):
            fit_params_out.append(NESTstats['marginals'][i]['median'])
            fit_params_out_std.append(NESTstats['marginals'][i]['sigma'])

        T,P,X = self.atmosphere.TP_profile(fit_params=fit_params_out)
        T_std = fit_params_out_std[self.fit_index[0]:self.fit_index[1]]
        X_std = fit_params_out_std[:self.fit_index[0]]

        return [(np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std))]

