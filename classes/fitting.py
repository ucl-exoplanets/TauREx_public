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

#conversion constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg


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
        #self.fit_params_init = self.fit_params
        #self.fit_index = self.atmosphere.fit_index


        # set priors, starting values, and fitting parameters
        self.set_fitting_data()

        #initialising output tags
        self.DOWNHILL = False
        self.MCMC = False
        self.NEST = False

        #DEBUG COUNTER
        self.db_count = 0

    def set_fitting_data(self):

        # todo this entire function can be removed. Feed priors and starting values directly in build_fit_params

        # set priors
        self.fit_priors = {
            'temp': (self.forwardmodel.atmosphere.planet_temp - self.params.fit_T_low,
                     self.forwardmodel.atmosphere.planet_temp + self.params.fit_T_up),
            # in AMU
            'mu': (self.params.fit_mu_low, self.params.fit_mu_up),
            # in RJUP
            'radius': (self.forwardmodel.atmosphere.planet_radius/RJUP - self.params.fit_radius_low,
                       self.forwardmodel.atmosphere.planet_radius/RJUP + self.params.fit_radius_up),
            # in Pascal
            'P0': (self.params.fit_P0_low, self.params.fit_P0_up),
            # 'X_absorbing': [], # todo maybe useless
            # 'X_inactive': [], # todo maybe useless
        }

        # # todo maybe useless
        # for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
        #     self.fit_priors['X_absorbing'].append((-np.inf,np.inf)) # working with centered-log-ratio transormation
        # for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
        #     self.fit_priors['X_inactive'].append((-np.inf,np.inf)) # working with centered-log-ratio transormation

        # set fitting parameters starting values
        self.fit_values = {
            'temp': self.forwardmodel.atmosphere.planet_temp,
            'mu': self.forwardmodel.atmosphere.planet_mu/AMU,
            'radius': self.forwardmodel.atmosphere.planet_radius/RJUP,
            'P0': np.mean(self.fit_priors['P0']),
            # 'X_absorbing': [], # todo maybe useless
            # 'X_inactive': [], # todo maybe useless
        }

        # todo maybe useless :
        # for idx, gasname in enumerate(self.inactive_gases):
        #     self.fit_values['X_absorbing'].append(0)
        # for n in range(len(self.inactive_gases)):
        #     self.fit_values['X_inactive'].append(0)

        # set fit_params list (input to fitting routines)
        self.build_fit_params()

    def build_fit_params(self):

        # build the fit_param list, input to all the fitting routines.
        # build the fit_bounds list, input to scipy.minimize and ?
        self.fit_params = []
        self.fit_params_names = []
        self.fit_bounds = []

        # include T, P0, r, mu only if we want to fit for them. Otherwise self.fit_values are assumed
        if not self.params.fit_fix_temp:
            self.fit_params.append(self.fit_values['temp'])
            self.fit_params_names.append('Temperature')
            self.fit_bounds.append(self.fit_priors['temp'])

        if not self.params.fit_couple_mu:
            if not self.params.fit_fix_mu:
                self.fit_params.append(self.fit_values['mu'])
                self.fit_params_names.append('mu')
                self.fit_bounds.append(self.fit_priors['mu'])

        if not self.params.fit_fix_radius:
            self.fit_params.append(self.fit_values['radius'])
            self.fit_params_names.append('Radius')
            self.fit_bounds.append(self.fit_priors['radius'])

        if not self.params.fit_fix_P0:
            self.fit_params.append(np.log(self.fit_values['P0']))
            self.fit_params_names.append('P0')
            self.fit_bounds.append(tuple(np.log(self.fit_priors['P0'])))

        # todo actually, we should set clr(x) = 0 for all ratios
        #
        # reparametrize the mixing ratios of absorbers and inactive gases using the centered-log-ratio transformation

        # get the geometrical mean
        sumlog = 0
        n = len(self.forwardmodel.atmosphere.absorbing_gases) + len(self.forwardmodel.atmosphere.inactive_gases)
        for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
            sumlog +=  np.log(self.forwardmodel.atmosphere.absorbing_gases_X[idx])

        for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
            sumlog += np.log(self.forwardmodel.atmosphere.inactive_gases_X[idx])
        mean = np.exp((1./n)*sumlog)

        # create array of clr(X + inactive_X)
        clr = []
        for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
            clr.append(np.log(self.forwardmodel.atmosphere.absorbing_gases_X[idx]/mean))
        for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
            clr.append(np.log(self.forwardmodel.atmosphere.inactive_gases_X[idx]/mean))

        # append clr(X) to the fitting value list
        self.fit_params += clr

        for val in clr:
            self.fit_bounds.append((-100, 100))

        # append molecule names to fit_params_names
        for gasname in self.forwardmodel.atmosphere.absorbing_gases:
            self.fit_params_names.append(gasname)
        for gasname in self.forwardmodel.atmosphere.inactive_gases:
            self.fit_params_names.append(gasname)

        self.fit_nparams = len(self.fit_params)


    def update_atmospheric_parameters(self, fit_params):

        # update atmospheric parameters with new values from fit_params

        count = 0

        # fitting for temperature
        if not self.params.fit_fix_temp:
            self.forwardmodel.atmosphere.planet_temp = fit_params[count]
            count += 1

        if not self.params.fit_couple_mu:
            # fitting for mean molecular weight mu
            if not self.params.fit_fix_mu:
                self.forwardmodel.atmosphere.planet_mu = fit_params[count]*AMU
                count += 1
        else:
            # coupling mu to atmospheric composition
            # planet_mu is a derived output, we're not directly fitting for mu...
            self.forwardmodel.atmosphere.planet_mu = self.forwardmodel.atmosphere.get_coupled_planet_mu()/AMU

        # fitting for the radius
        if not self.params.fit_fix_radius:
            self.forwardmodel.atmosphere.planet_radius = fit_params[count]*RJUP
            count += 1

        # we need to update surface gravity and scale height if temp, mu, or radius have changed
        self.forwardmodel.atmosphere.planet_grav = self.forwardmodel.atmosphere.get_surface_gravity()
        self.forwardmodel.atmosphere.scaleheight = self.forwardmodel.atmosphere.get_scaleheight()

        # fitting for surface pressure
        if not self.params.fit_fix_P0:
            self.forwardmodel.atmosphere.max_pressure = np.exp(fit_params[count])
            count += 1

        # update TP profile
        self.forwardmodel.atmosphere.pta = self.forwardmodel.atmosphere.setup_pta_grid()
        self.forwardmodel.atmosphere.P = self.forwardmodel.atmosphere.pta[:,0] # pressure array
        self.forwardmodel.atmosphere.P_bar = self.forwardmodel.atmosphere.P * 1.0e-5 #convert pressure from Pa to bar
        self.forwardmodel.atmosphere.T = self.forwardmodel.atmosphere.pta[:,1] # temperature array
        self.forwardmodel.atmosphere.z = self.forwardmodel.atmosphere.pta[:,2] # altitude array
        self.forwardmodel.atmosphere.rho = self.forwardmodel.atmosphere.get_rho() # update density

        # @todo mixing ratios are always fitted. Maybe allow to fix, either absorbing and/or inactive gases

        # fitting for mixing ratios. Need to convert clr(X) to X (centerd-log-ratio inverse transformation)

        # get the geometric mean of all mixing ratios
        ngases = len(self.forwardmodel.atmosphere.absorbing_gases) + len(self.forwardmodel.atmosphere.inactive_gases)
        clr_inv = []
        for i in range(ngases):
            clr_inv.append(np.exp(fit_params[count+i]))
        clr_inv = np.asarray(clr_inv)/np.sum(clr_inv)

        count2 = 0
        for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
            self.forwardmodel.atmosphere.absorbing_gases_X[idx] = clr_inv[count2]
            count2 += 1
        for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
            self.forwardmodel.atmosphere.inactive_gases_X[idx] = clr_inv[count2]
            count2 += 1

        # @todo careful, next line won't work if preselector is running. Fix preselector!
        self.forwardmodel.atmosphere.X = self.forwardmodel.atmosphere.set_mixing_ratios() # mixing ratios are read from parameter file or set to 1e-4 if preselector = True

        # recalculate Rayleigh scattering
        if self.params.in_include_Rayleigh:
            self.forwardmodel.Rsig = self.forwardmodel.get_Rsig()
        else:
            self.forwardmodel.Rsig = zeros((self.forwardmodel.nlambda))


    def chisq_trans(self, fit_params, data, datastd):

        # chisquare minimisation bit
        # only working for transmission

        # update atmospheric parameters with fit_params values
        self.update_atmospheric_parameters(fit_params)

        # get forward model
        model = self.forwardmodel.model()

        #binning internal model
        model_binned = [model[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid)]
        #
        # clf()
        # errorbar(self.data.spectrum[:,0],self.data.spectrum[:,1], yerr=self.data.spectrum[:,2])
        # plot(self.data.spectrum[:,0], model_binned)
        # draw()
        # pause(0.00001)

        res = (data - model_binned) / datastd

        # print '%i  |   %.1f, %.6f, %.4f, %.4f' % (sqrt(sum(res**2)), self.forwardmodel.atmosphere.planet_temp, \
        #     self.forwardmodel.atmosphere.planet_mu, \
        #     self.forwardmodel.atmosphere.planet_radius/RJUP, \
        #     self.forwardmodel.atmosphere.max_pressure/1.e5), \
        #     self.forwardmodel.atmosphere.absorbing_gases_X, \
        #     self.forwardmodel.atmosphere.inactive_gases_X

        return sum(res**2)

###############################################################################
#simplex downhill algorithm

    def downhill_fit(self):

        ion()

        logging.info('Fit data using %s minimisation' % self.params.downhill_type)

        # @todo here we could use different starting points for different threads?

        data = self.data.spectrum[:,1] # observed data
        datastd = self.data.spectrum[:,2] # data error

        fit_params = self.build_fit_params()

        fit_output = minimize(fun=self.chisq_trans,
                              x0=self.fit_params,
                              args=(data,datastd),
                              method=self.params.downhill_type,
                              bounds=(self.fit_bounds))

        logging.info('Saving the downhill minimization results')

        self.update_atmospheric_parameters(fit_output['x'])

        self.DOWNHILL = True
        self.DOWNHILL_T_mean = self.forwardmodel.atmosphere.planet_temp
        self.DOWNHILL_X_mean = self.forwardmodel.atmosphere.X
        self.DOWNHILL_fit_output = fit_output['x']

###############################################################################    
#Markov Chain Monte Carlo algorithm


    # noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment
    # @profile #line-by-line profiling decorator
    def mcmc_fit(self):
    # adaptive Markov Chain Monte Carlo

        ion()


        logging.info('Start MCMC fit')

        if self.DOWNHILL:
            fit_params_init = self.DOWNHILL_fit_output
        else:
            fit_params_init = self.fit_params

        data = self.data.spectrum[:,1] #observed data
        datastd = self.data.spectrum[:,2] #data error
        
        # setting prior distributions
        # master thread (MPIrank =0) will start from ideal solution (from downhill fitting)
        # slave threads (MPIrank != 0) will start from randomised points.

        priors = empty(size(fit_params_init), dtype=object)

        #setting up main thread. Use downhill FIT as starting points
        if self.MPIrank == 0:
            for i in range(self.fit_nparams):
                    priors[i] = pymc.Uniform('PFIT_%i' % (i), self.fit_bounds[i][0], self.fit_bounds[i][1],
                                             value=fit_params_init[i])  # uniform prior

        else:
            #setting up other threads (if exist). Their initial starting positions will be randomly perturbed
            for i in range(self.fit_nparams):
                param_range = (selffit_bounds[i][1] - self.fit_bounds[i][0]) / 5.0 #range of parameter over which to perturb starting position
                param_mean  = np.mean(self.fit_bounds[i])
                param_rand  = random.uniform(low=param_mean-param_range,high=param_mean+param_range) #random parameter start
                priors[i] = pymc.Uniform('PFIT_%i' % (i), self.fit_bounds[i][0],self.fit_bounds[i][1],value=param_rand)  # uniform prior

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

        T,P,X = self.forwardmodel.atmosphere.TP_profile(fit_params=fit_params_out)

        T_std = fit_params_out_std[self.fit_index[0]:self.fit_index[1]]
        X_std = fit_params_out_std[:self.fit_index[0]]

        return np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std)


############################################################################### 
#Nested Sampling algorithm

    def multinest_fit(self, resume=None):
        #multinest fitting routine (wrapper around PyMultiNest)


        ion()


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
            fit_params_container = asarray([cube[i] for i in xrange(len(self.fit_params))])
            chi_t = self.chisq_trans(fit_params_container, data, datastd)
            llterms = (-ndim/2.0)*np.log(2.*pi*datastd_mean**2) - 0.5*chi_t
            return llterms
    
        def multinest_uniform_prior(cube, ndim, nparams):
            #prior distributions called by multinest. Implements a uniform prior
            #converting parameters from normalised grid to uniform prior
            for i in xrange(len(self.fit_params)):
                cube[i] = (cube[i] * (self.fit_bounds[i][1]-self.fit_bounds[i][0])) + self.fit_bounds[i][0]

        data = self.data.spectrum[:,1] #observed data
        datastd = self.data.spectrum[:,2] #data error
        datastd_mean = mean(datastd)
        
        parameters = [str(i) for i in xrange(len(self.fit_params))]
        ndim = len(self.fit_params)

        #progress = pymultinest.ProgressPlotter(n_params = len(self.fit_params)); progress.start()
        #threading.Timer(60, show, ["chains/1-phys_live.points.pdf"]).start() # delayed opening
        logging.info('Multinest output dir is %s' % self.dir_multinest)

        pymultinest.run(LogLikelihood=multinest_loglike,
                        Prior=multinest_uniform_prior,
                        n_dims=len(self.fit_params),
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
            NESTout = pymultinest.Analyzer(n_params=len(self.fit_params),
                                           outputfiles_basename=os.path.join(self.dir_multinest, '1-'))


            # this returns a list of output solutions: [[fit_params_mode1], [fit_params_mode2], ... ]
            self.multinest_result = self.collate_multinest_result(NESTout)


            #saving arrays to object
            self.NEST = True
            self.NEST_modes = len(self.multinest_result) # number of modes detected

            # self.NEST_T_mean = [solution[0] for solution in multinest_result]
            # self.NEST_T_std = [solution[1] for solution in multinest_result]
            # self.NEST_X_mean = [solution[2] for solution in multinest_result]
            # self.NEST_X_std = [solution[3] for solution in multinest_result]

            self.NEST_stats = NESTout.get_stats()
            self.NEST_FITDATA = NESTout

    def collate_multinest_result(self, NESTout):

        # Multinest in multimdes can return two or more sets of solutions.
        # if only one solution is found, it returns a list with a tuple with 4 elements: [(T, T_std, X, X_std)]
        # if N (>1) solutions are found, it returns a list with N tuples with 4 elements each
        #    e.g. [(T1, T_std1, X1, X_std1), (T2 T_std2, X2, X_std2), ... ]


        logging.info('Unpacking the MULTINEST results')

        # todo: standard deviations are not stored...
        # todo: oh wow I did a mess! Improve managment of outputs

        NESTstats = NESTout.get_stats()
        out_list = []

        if self.params.nest_multimodes:
            if len(NESTstats['modes']) > 1:
                # more than one mode found
                for n in range(len(NESTstats['modes'])):
                    fit_params_out = []
                    fit_params_out_std = []
                    for i in range(len(self.fit_params)):
                        fit_params_out.append(NESTstats['modes'][n]['maximum a posterior'][i])
                        #fit_params_out_std.append(NESTstats['modes'][n]['sigma'][i])
                    out_list.append(fit_params_out)
        else:
            # only one set of solutions found
            fit_params_out = []
            fit_params_out_std = []
            for i in range(len(self.fit_params)):
                fit_params_out.append(NESTstats['marginals'][i]['median'])
                #fit_params_out_std.append(NESTstats['marginals'][i]['sigma'])
            out_list.append(fit_params_out)

        return out_list
