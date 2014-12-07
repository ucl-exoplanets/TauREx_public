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
#    def __init__(self, params, data, profile, rad_model=None):
    def __init__(self, profile, params=None, data=None, rad_model=None):

        # @todo here we should actually have the transmission or emission objects as input. Then get rid of self.set_model()

        logging.info('Initialise object fitting')

        if params:
            self.params = params
        else:
            self.params = profile.params #get params object from profile

        if data:
            self.data = data
        else:
            self.data = profile.data    # get data object from profile

        self.profile     = profile

        # this is equivalent to type(instance)
        self.__ID__      = 'fitting' 
        
        # MPI support
        if MPIimport:
            self.MPIrank     = MPI.COMM_WORLD.Get_rank()
            self.MPIsize     = MPI.COMM_WORLD.Get_size()
        else:
            self.MPIrank     = 0
            self.MPIsize     = 0

        # set some folder names
        self.dir_chains = os.path.join(self.params.out_path, 'chains')
        self.dir_mcmc = os.path.join(self.dir_chains, 'MCMC')
        self.dir_mutlinest = os.path.join(self.dir_chains, 'multinest')

        # create some folders.
        if self.MPIrank == 0:
            folders = [self.params.out_path, self.dir_chains, self.dir_mcmc, self.dir_mutlinest]
            for f in folders:
                if not os.path.isdir(f):
                    logging.info('Create folder %s' % f)
                    os.mkdir(f)

        #loading transmission/emission model
        self.set_model(rad_model)

        #loading data and parameters
        self.observation = self.data.spectrum
        self.nlayers     = self.params.tp_atm_levels # @todo nlayers should be stored in the tp_profile object
        self.ngas        = self.data.ngas
        self.nspecbingrid= len(self.data.spec_bin_grid) # @todo move to data object

        #self.n_params    = int(self.ngas  + 1) #+1 for only one temperature so far
        #self.n_params    = 2 #restricting to one temperature and one column density parameter at the moment

        #defining prior bounds  @todo remove? These are now in self.bounds
        self.X_up        = self.params.fit_X_up
        self.X_low       = self.params.fit_X_low
        self.T_up        = self.params.planet_temp + self.params.fit_T_up
        self.T_low       = self.params.planet_temp - self.params.fit_T_low

        #getting parameter grid with initial estimates
        self.PINIT       = self.profile.PARAMS  # @todo find better var name (confusion with params!)
        self.Pindex      = self.profile.TPindex
        self.Pshape      = shape(self.PINIT)
        self.n_params    = len(self.PINIT) # @todo find better var name (confusion with params!)
        
        #getting prior bounds for downhill algorithm
        self.bounds      = self.profile.bounds

        #initialising output tags
        self.DOWNHILL = False
        self.MCMC     = False
        self.NEST     = False

        #DEBUG COUNTER
        self.db_count = 0

    # @profile #line-by-line profiling decorator
    def chisq_trans(self,PFIT,DATA,DATASTD):

        #chisquare minimisation bit

        T,P,X = self.profile.TP_profile(PARAMS=PFIT)  # calculating TP profile

        rho = self.profile.get_rho(T=T, P=P)          # calculating densities

        #the temperature parameter should work out of the box but check for transmission again
        MODEL = self.model(rho=rho, X=X, temperature=T) # model() is in the base class!

#         MODEL = self.transmod.cpath_integral(rho=rho,X=X,temperature=1400)

        #binning internal model
        MODEL_binned = [MODEL[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.nspecbingrid)]


        #         MODEL_interp = np.interp(self.data.wavegrid,self.data.specgrid,MODEL)
#         print PFIT
#         print T
#         print P
#         ion()
#         clf()
#         figure(1)
#         plot(T,np.log(P))
#         draw()

        #print PFIT, T, P
#        plt.clf()
#        plt.plot(DATA,'g')
#        plt.plot(MODEL_binned)
#        plt.draw()
#        plt.pause(0.0001)

        # MODEL_interp = DATA+np.random.random(np.shape(MODEL_interp))
        
        res = (DATA-MODEL_binned) / DATASTD
#         print sum(res**2)
        return sum(res**2)
        
############################################################################### 
#simplex downhill algorithm

    # @profile #line-by-line profiling decorator
    def collate_downhill_results(self,PFIT):

        logging.info('Unpacking the downhill minimization results')

        T,P,X = self.profile.TP_profile(PARAMS=PFIT)

        return T,X 
    
    # @profile #line-by-line profiling decorator
    def downhill_fit(self):

        logging.info('Fits data using simplex downhill minimisation')

  #      plt.ion()
  #      plt.show()

        # @todo here we could use different starting points for different threads?
        PINIT = self.PINIT  # initial temperatures and abundances
        DATA = self.observation[:,1] # observed data
        DATASTD = self.observation[:,2] # data error
        
        # PFIT, err, out3, out4, out5 = fmin(self.chisq_trans, PINIT, args=(DATA,DATASTD), xtol=1e-5, ftol=1e-5,maxiter=1e6,
        #                                    disp=1, full_output=1)

        PFIT = minimize(fun=self.chisq_trans,
                        x0=PINIT,
                        args=(DATA,DATASTD),
                        method=self.params.downhill_type,
                        bounds=(self.bounds))
#         PFIT = minimize(self.chisq_trans,PINIT,args=(DATA,DATASTD),method='Nelder-Mead',bounds=(self.bounds))
#         PFIT = fmin(self.chisq_trans,PINIT,args=(DATA,DATASTD),maxfun=10)

        Tout_mean, Xout_mean = self.collate_downhill_results(PFIT['x'])

        self.DOWNHILL = True
        self.DOWNHILL_T_mean = Tout_mean
        self.DOWNHILL_X_mean = Xout_mean
        self.DOWNHILL_PFIT   = PFIT['x']

###############################################################################    
#Markov Chain Monte Carlo algorithm

    # @profile #line-by-line profiling decorator
    def collate_mcmc_result(self,MCMCout):

        logging.info('Unpacking the MCMC results')

        MCMCstats = MCMCout.stats()
        
        PFIT     = []
        PFIT_std = []
        for i in range(len(self.PINIT)):
            PFIT.append(MCMCstats['PFIT_%i' % i]['mean'])
            PFIT_std.append(MCMCstats['PFIT_%i' % i]['standard deviation'])
            
        T,P,X = self.profile.TP_profile(PARAMS=PFIT)

        T_std = PFIT_std[self.Pindex[0]:self.Pindex[1]]
        X_std = PFIT_std[:self.Pindex[0]]

        return np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std)


    # noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment
    # @profile #line-by-line profiling decorator
    def mcmc_fit(self):
    # adaptive Markov Chain Monte Carlo

        logging.info('Start MCMC fit')

        if self.DOWNHILL:
            PINIT = self.DOWNHILL_PFIT
        else:
            PINIT = self.PINIT
        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        
        
        # setting prior distributions
        # master thread (MPIrank =0) will start from ideal solution (from downhill fitting)
        # slave threads (MPIrank != 0) will start from randomised points.

        priors = empty(size(PINIT),dtype=object)

        #setting up main thread. Use downhill FIT as starting points
        if self.MPIrank == 0:
            for i in range(self.n_params):
                    priors[i] = pymc.Uniform('PFIT_%i' % (i), self.bounds[i][0],self.bounds[i][1],value=PINIT[i])  # uniform prior

        #setting up other threads (if exist). Their initial starting positions will be randomly perturbed
        else:
            for i in range(self.n_params):
                P_range = (self.bounds[i][1] - self.bounds[i][0]) / 5.0 #range of parameter over which to perturb starting position
                P_mean  = np.mean(self.bounds[i])
                P_rand  = random.uniform(low=P_mean-P_range,high=P_mean+P_range) #random parameter start
#                print self.bounds[i][0], self.bounds[i][1], P_mean, P_range, P_rand
                priors[i] = pymc.Uniform('PFIT_%i' % (i), self.bounds[i][0],self.bounds[i][1],value=P_rand)  # uniform prior

        #setting up data error prior if specified
        if self.params.mcmc_update_std:
            std_dev = pymc.Uniform('std_dev',0.0,2.0*max(DATASTD),value=DATASTD,size=len(DATASTD)) #uniform prior on data standard deviation
        else:
            std_dev = pymc.Uniform('std_dev',0.0,2.0*max(DATASTD),value=DATASTD,observed=True,size=len(DATASTD))
              

#         @pymc.stochastic(observed=True, plot=False)
#         def mcmc_loglikelihood(value=PINIT, PFIT=priors, DATASTD=precision, DATA=DATA):
        def mcmc_loglikelihood(value, PFIT, DATASTD, DATA):
            # log-likelihood function. Needs to be initialised directly since CYTHON does not like PYMC decorators
            # @todo need to cq ast from numpy object array to float array. Slow but hstack is  slower.

            Pcontainer = np.zeros((self.n_params))
            for i in range(self.n_params):
                Pcontainer[i] = PFIT[i]

            # @todo Pcontainer should be equal to PFIT? I think so...

            chi_t = self.chisq_trans(Pcontainer, DATA, DATASTD) #calculate chisq

            llterms =   (-len(DATA)/2.0)*np.log(pi) - np.log(np.mean(DATASTD)) - 0.5* chi_t

            return llterms

        mcmc_logp = pymc.Stochastic( logp = mcmc_loglikelihood,
                doc = 'The switchpoint for mcmc loglikelihood.',
                name = 'switchpoint',
                parents = {'PFIT': priors, 'DATASTD': std_dev, 'DATA': DATA},
#                 random = switchpoint_rand,
                trace = True,
                value = PINIT,
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
        self.MCMC        = True
        self.MCMC_T_mean = Tout_mean
        self.MCMC_T_std  = Tout_std
        self.MCMC_X_mean = Xout_mean
        self.MCMC_X_std  = Xout_std
        self.MCMC_STATS  = R.stats()
        self.MCMC_FITDATA= R

        # each MCMC chain is run on a separate thread. Save all outputs
        logging.info('Store the MCMC results')
        self.MCMC_OUT = {}
        self.MCMC_OUT[self.MPIrank] = {}
        self.MCMC_OUT[self.MPIrank]['FITDATA'] = R
        self.MCMC_OUT[self.MPIrank]['STATS']   = R.stats()
        self.MCMC_OUT[self.MPIrank]['T_mean']  = Tout_mean
        self.MCMC_OUT[self.MPIrank]['T_std']   = Tout_std
        self.MCMC_OUT[self.MPIrank]['X_mean']  = Xout_mean
        self.MCMC_OUT[self.MPIrank]['X_std']   = Xout_std

        #saving MCMC results to file as pickle. Do it later.
        #with gzip.GzipFile(self.params.out_path+'MCMC_results.pkl.zip','wb') as outhandle:
        #     pickle.dump(MCMC_OUT,outhandle)


############################################################################### 
#Nested Sampling algorithm

    def collate_multinest_result(self,NESTout):

        logging.info('Unpacking the MULTINEST results')

        NESTstats = NESTout.get_stats()

        PFIT     = []
        PFIT_std = []
      
        for i in range(self.n_params):
            if self.params.nest_multimodes:
                PFIT.append(NESTstats['modes'][0]['maximum a posterior'][i])
                PFIT_std.append(NESTstats['modes'][0]['sigma'][i])
            else:
                PFIT.append(NESTstats['marginals'][i]['median'])
                PFIT_std.append(NESTstats['marginals'][i]['sigma'])

        T,P,X = self.profile.TP_profile(PARAMS=PFIT)
        
        T_std = PFIT_std[self.Pindex[0]:self.Pindex[1]]
        X_std = PFIT_std[:self.Pindex[0]]
        # Xout_mean = Xout_mean.reshape(self.ngas,self.nlayers)
        # Xout_std = Xout_std.reshape(self.ngas,self.nlayers)

        return np.asarray(T), np.asarray(T_std), np.asarray(X), np.asarray(X_std)


    
    def multinest_fit(self,resume=None):
        #multinest fitting routine (wrapper around PyMultiNest)

        logging.info('Start MULTINEST fit')

        if resume is None:
            resume = self.params.nest_resume

        if resume:
            filename = os.path.join(self.dir_mutlinest, '1-live.points')
            if os.path.isfile(filename):
                livepoints = sum(1 for line in open(filename))
                if livepoints != self.params.nest_nlive:
                    logging.warning('Cannot resume previous MULTINEST run, the number of live points has changed')
                    logging.warning('Removing previous MULTINEST chains')
                    for file in os.listdir(self.dir_mutlinest):
                        file_path = os.path.join(self.dir_mutlinest, file)
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

        def multinest_loglike(cube, ndim,nparams):
            #log-likelihood function called by multinest
            PFIT = [cube[i] for i in xrange(self.n_params)]
            PFIT = asarray(PFIT)
            chi_t = self.chisq_trans(PFIT,DATA,DATASTD)
            llterms =   (-ndim/2.0)*np.log(2.*pi*DATASTDmean**2) - 0.5*chi_t
            return llterms
    
        def multinest_uniform_prior(cube,ndim,nparams):
            #prior distributions called by multinest. Implements a uniform prior

            #converting parameters from normalised grid to uniform prior
            for i in xrange(self.n_params):
                cube[i] = (cube[i] * (self.bounds[i][1]-self.bounds[i][0])) + self.bounds[i][0]
#            cube[1] = cube[1] * (self.X_up-self.X_low) + self.X_low

        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        DATASTDmean = mean(DATASTD)
        
        parameters = [str(i) for i in xrange(self.n_params)]
        n_params = self.n_params
        ndim = n_params

        #progress = pymultinest.ProgressPlotter(n_params = n_params); progress.start()
        #threading.Timer(60, show, ["chains/1-phys_live.points.pdf"]).start() # delayed opening
        print 'Livepoints %i ' % self.params.nest_nlive
        pymultinest.run(LogLikelihood=multinest_loglike,
                        Prior=multinest_uniform_prior,
                        n_dims=self.n_params,
                        multimodal=self.params.nest_multimodes,
                        max_modes=self.params.nest_max_modes,
                        outputfiles_basename=os.path.join(self.dir_mutlinest, '1-'),
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
            OUT = pymultinest.Analyzer(n_params=self.n_params)

            Tout_mean, Tout_std, Xout_mean, Xout_std = self.collate_multinest_result(OUT)

            logging.info('Store the MULTINEST results')

            #saving arrays to object
            self.NEST        = True
            self.NEST_T_mean = Tout_mean
            self.NEST_T_std  = Tout_std
            self.NEST_X_mean = Xout_mean
            self.NEST_X_std  = Xout_std
            self.NEST_STATS  = OUT.get_stats()
            self.NEST_FITDATA= OUT