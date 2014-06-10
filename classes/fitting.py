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
import numpy, pylab,os,sys,math,pymc,warnings,threading, subprocess,gzip,pickle
from numpy import *
from pylab import *
from StringIO import StringIO
from scipy.interpolate import interp1d
from scipy.optimize import fmin
from scipy.optimize import minimize

try: 
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

if not os.path.exists("chains"): os.mkdir("chains")

try:
    from mpi4py import MPI
except ImportError:
    pass




class fitting(object):
    def __init__(self,params,data,profile,transmod):
        
        self.MPIrank     = MPI.COMM_WORLD.Get_rank()
        
        self.params      = params
        self.dataob      = data
        self.transmod    = transmod
        self.profileob   = profile
        self.observation = data.spectrum
        self.nlayers     = params.tp_atm_levels
        self.ngas        = data.ngas
        self.n_params    = int(self.ngas  + 1) #+1 for only one temperature so far
        
        # self.n_params    = 2 #restricting to one temperature and one column density parameter at the moment

        #defining prior bounds
        self.X_up        = self.params.fit_X_up
        self.X_low       = self.params.fit_X_low
        self.T_up        = self.params.planet_temp + self.params.fit_T_up
        self.T_low       = self.params.planet_temp - self.params.fit_T_low

        #setting bounds for downhill algorithm
        self.bounds = []
        self.bounds.append((self.T_low,self.T_up))
        for i in range(self.ngas):
            self.bounds.append((self.X_low,self.X_up))
            


        #checking if number of free parameters for X = number of gas species
#         if len(self.params.fit_param_free_X) != (self.params.tp_atm_levels*self.params.tp_num_gas):
#             warnings.warn('Number of free X parameters does not equal number of gas species and atmospheric levels')
#             self.params.fit_param_free_X = self.params.tp_atm_levels*self.params.tp_num_gas
        
        #setting up initial parameter array
        # self.PINIT = zeros((self.ngas+1,self.nlayers))
        # self.PINIT[0,:] = self.profileob.T
        # self.PINIT[1:,:] = self.profileob.X

        self.PINIT     = zeros((self.ngas+1))
        self.PINIT[0]  = self.profileob.T[0]
        # self.PINIT[1:] = float((self.params.fit_X_up-self.params.fit_X_low)/2 + self.params.fit_X_low)
        self.PINIT[1:] = 1e-5


        # for i in range(1,self.ngas):
        #     self.PINIT[i] = self.profileob.X[i,0]

        self.Pshape = shape(self.PINIT)
        
        # self.PINIT2 = self.PINIT.flatten(order='F')
        # self.PINIT3 = self.PINIT2[(self.Pshape[0]-1):] #restricting temperature to one free parameter only
        self.PINIT3 = self.PINIT

       
        #initialising output tags
        self.DOWNHILL = False
        self.MCMC     = False
        self.NEST     = False



    # @profile #line-by-line profiling decorator
    def chisq_trans(self,PFIT,DATA,DATASTD):
        #chisquare minimisation bit

        rho = self.profileob.get_rho(T=PFIT[0])
#         print PFIT[0]
#         X   = PFIT[1:].reshape(self.Pshape[0],self.Pshape[1]-1)

        X   = zeros((self.ngas,self.nlayers))
        for i in range(self.ngas):
            X[i,:] += PFIT[i+1]


        MODEL = self.transmod.cpath_integral(rho=rho,X=X,temperature=PFIT[0])
        MODEL_interp = np.interp(self.dataob.wavegrid,self.dataob.specgrid,MODEL)

#         ion()
#         clf()
#         figure(100)
#         plot(MODEL)
#         draw()
#         clf()
#         figure(101)
#         plot(DATA,'g')
#         plot(MODEL_interp)
#         draw()

        # MODEL_interp = DATA+np.random.random(np.shape(MODEL_interp))

        res = (DATA-MODEL_interp) / DATASTD
        return sum(res**2)
        
############################################################################### 
#simplex downhill algorithm

    # @profile #line-by-line profiling decorator
    def collate_downhill_results(self,PFIT):
        #function unpacking the downhill minimization results
        Xout_mean = zeros((self.ngas,self.nlayers))
        # Tout_mean = zeros((self.nlayers*self.ngas))

        print PFIT
        Tout_mean= PFIT[0]

        if len(PFIT)-1 < self.nlayers*self.ngas:
            # Xout_mean[:] += PFIT[1:]
            for i in range(len(PFIT)-1):
                Xout_mean[i,:] += PFIT[i+1]
        else:
            Xout_mean[:] = PFIT[1:]

        return Tout_mean,Xout_mean

    # @profile #line-by-line profiling decorator
    def downhill_fit(self):
    # fits data using simplex downhill minimisation

        PINIT = self.PINIT3  # initial temperatures and abundances
        
        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        
        # PFIT, err, out3, out4, out5 = fmin(self.chisq_trans, PINIT, args=(DATA,DATASTD), xtol=1e-5, ftol=1e-5,maxiter=1e6,
        #                                    disp=1, full_output=1)

        PFIT = minimize(self.chisq_trans,PINIT,args=(DATA,DATASTD),method='L-BFGS-B',bounds=(self.bounds))


        Tout_mean, Xout_mean = self.collate_downhill_results(PFIT['x'])
        self.DOWNHILL = True
        self.DOWNHILL_T_mean = Tout_mean
        self.DOWNHILL_X_mean = Xout_mean
        self.DOWNHILL_PFIT   = PFIT['x']


    
###############################################################################    
#Markov Chain Monte Carlo algorithm

    # @profile #line-by-line profiling decorator
    def collate_mcmc_result(self,MCMCout):
        #function unpacking the MCMC results
        
        MCMCstats = MCMCout.stats()
        Xout_mean = zeros((self.ngas,self.nlayers))
        Xout_std  = zeros((self.ngas,self.nlayers))

        #needs to be rewritten to make assignment more dynamic for variable TP profiles.Ideas?
        Tout_mean = MCMCstats['temp']['mean']
        Tout_std  = MCMCstats['temp']['standard deviation']

        # print MCMCstats.keys()
        # if len(MCMCstats.keys())-1 < self.nlayers*self.ngas:
        #     Xout_mean[:] += MCMCstats['mixing_%i' % 0]['mean']
        #     Xout_std[:]  += MCMCstats['mixing_%i' % 0]['standard deviation']
        #     for i in range(len(MCMCstats.keys())-1):
        #         Xout_mean[i] = MCMCstats['mixing_%i' % i]['mean']
        #         Xout_std[i]  = MCMCstats['mixing_%i' % i]['standard deviation']
        # else:
        #     for i in range(self.nlayers*self.ngas):
        #         Xout_mean[i] = MCMCstats['mixing_%i' % i]['mean']
        #         Xout_std[i]  = MCMCstats['mixing_%i' % i]['standard deviation']

        for i in range(self.ngas):
            Xout_mean[i,:] = MCMCstats['mixing_%i' % i]['mean']
            Xout_std[i,:]  = MCMCstats['mixing_%i' % i]['standard deviation']

        # Xout_mean = Xout_mean.reshape(self.ngas,self.nlayers)
        # Xout_std = Xout_std.reshape(self.ngas,self.nlayers)

        return Tout_mean, Tout_std, Xout_mean, Xout_std


    # noinspection PyNoneFunctionAssignment,PyNoneFunctionAssignment
    # @profile #line-by-line profiling decorator
    def mcmc_fit(self):
    # adaptive Markov Chain Monte Carlo

        if self.DOWNHILL:
            PINIT = self.DOWNHILL_PFIT
        else:
            PINIT = self.PINIT3
        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        
        
        # setting prior distributions
        priors = empty(size(PINIT),dtype=object)
        priors[0] = pymc.Uniform('temp',self.T_low,self.T_up,value=PINIT[0]) #uniform temperature prior
        for i in range(len(PINIT)-1):
            priors[i+1] = pymc.Uniform('mixing_%i' % (i), self.X_low,self.X_up,value=PINIT[i+1])  # uniform mixing ratio prior
        

        #setting up data error prior if specified
        if self.params.mcmc_update_std:
            std_dev = pymc.Uniform('std_dev',0.0,max(DATASTD),value=DATASTD,size=len(DATASTD)) #uniform prior on data standard deviation
        else:
            std_dev = pymc.Uniform('std_dev',0.0,max(DATASTD),value=DATASTD,observed=True,size=len(DATASTD))
        
        
        @pymc.deterministic(plot=False)
        def precision(std_dev=std_dev):
                return std_dev
    
        
        # log-likelihood function
        @pymc.stochastic(observed=True, plot=False)
        def mcmc_loglikelihood(value=PINIT, PFIT=priors, DATASTD=precision, DATA=DATA):
            
            chi_t = self.chisq_trans(PFIT,DATA,DATASTD)
            # llterms = len(DATA) * np.log(1.0 / sqrt(2.0 * pi)) - len(DATA) * np.log(DATASTD) - 0.5 * chi_t
            llterms =  - 0.5* chi_t
            return llterms

        # executing MCMC sampling
        R = pymc.MCMC((priors, mcmc_loglikelihood), verbose=1)  # build the model
        R.sample(iter=self.params.mcmc_iter, burn=self.params.mcmc_burn,
                 thin=self.params.mcmc_thin)              # populate and run it
        
        #coallating results into arrays
        Tout_mean, Tout_std, Xout_mean, Xout_std = self.collate_mcmc_result(R)

        print Xout_mean, Tout_mean
        #saving arrays to object
        self.MCMC        = True
        self.MCMC_T_mean = Tout_mean
        self.MCMC_T_std  = Tout_std
        self.MCMC_X_mean = Xout_mean
        self.MCMC_X_std  = Xout_std
        self.MCMC_STATS  = R.stats()
        self.MCMC_FITDATA= R
        
        MCMC_OUT = {}
        MCMC_OUT[self.MPIrank] = {}
        MCMC_OUT[self.MPIrank]['FITDATA'] = R
        MCMC_OUT[self.MPIrank]['STATS']   = R.stats()
        MCMC_OUT[self.MPIrank]['T_mean']  = Tout_mean
        MCMC_OUT[self.MPIrank]['T_std']   = Tout_std
        MCMC_OUT[self.MPIrank]['X_mean']  = Xout_mean
        MCMC_OUT[self.MPIrank]['X_std']   = Xout_std
        
        
        #saving to file as pickle 
        with gzip.GzipFile(self.params.out_path+'MCMC_results.pkl.zip','wb') as outhandle:
            pickle.dump(MCMC_OUT,outhandle)



    
############################################################################### 
#Nested Sampling algorithm

    def collate_multinest_result(self,NESTout):
        #function unpacking the MULTINEST results

        NESTstats = NESTout.get_stats()

        Xout_mean = zeros((self.ngas,self.nlayers))
        Xout_std  = zeros((self.ngas,self.nlayers))

        # print NESTstats['modes'][0]['maximum a posterior']
        Tout_mean = NESTstats['modes'][0]['maximum a posterior'][0]
        Tout_std  = NESTstats['modes'][0]['sigma'][0]

        # for i in range(1,len(NESTstats['modes'][0]['maximum a posterior'])):
        #     Xout_mean[i] = NESTstats['modes'][0]['maximum a posterior'][i]
        #     Xout_std[i]  = NESTstats['modes'][0]['sigma'][i]
        #
        # if len(NESTstats['modes'][0]['maximum a posterior'])-1 < self.nlayers*self.ngas:
        #     Xout_mean[:] += NESTstats['modes'][0]['maximum a posterior'][1]
        #     Xout_std[:]  += NESTstats['modes'][0]['sigma'][1]
        #     for i in range(len(NESTstats['modes'][0]['maximum a posterior'])-1):
        #         Xout_mean[i] = NESTstats['modes'][0]['maximum a posterior'][i]
        #         Xout_std[i]  = NESTstats['modes'][0]['sigma'][i]
        # else:
        #     for i in range(self.nlayers*self.ngas):
        #         Xout_mean[i] = NESTstats['modes'][0]['maximum a posterior'][i]
        #         Xout_std[i]  = NESTstats['modes'][0]['sigma'][i]

        for i in range(self.ngas):
            Xout_mean[i,:] = NESTstats['modes'][0]['maximum a posterior'][i+1]
            Xout_std[i,:]  = NESTstats['modes'][0]['sigma'][i+1]


        # Xout_mean = Xout_mean.reshape(self.ngas,self.nlayers)
        # Xout_std = Xout_std.reshape(self.ngas,self.nlayers)

        return Tout_mean, Tout_std, Xout_mean, Xout_std

    
    def multinest_fit(self,resume=None):
        #multinest fitting routine (wrapper around PyMultiNest)
        if resume is None:
            resume = self.params.nest_resume

        def show(filepath): 
            """ open the output (pdf) file for the user """
            if os.name == 'mac': subprocess.call(('open', filepath))
            elif os.name == 'nt': os.startfile(filepath)
        
        
        
        def multinest_loglike(cube, ndim,nparams):
        #log-likelihood function called by multinest
        
            PFIT = [cube[i] for i in range(self.n_params)]
            PFIT = asarray(PFIT)

            chi_t = self.chisq_trans(PFIT,DATA,DATASTD)
            llterms =   (-ndim/2.0)*np.log(pi) -np.log(DATASTDmean) -0.5* chi_t
#             llterms =    -0.5* chi_t
            return llterms
    
    
        def multinest_uniform_prior(cube,ndim,nparams):
            #prior distributions called by multinest
            #implements a uniform prior

            #converting temperatures from normalised grid to uniform prior
            cube[0] = cube[0]* (self.T_up-self.T_low)+self.T_low
    #         print cube[0]

            #converting mixing ratios from normalised grid to uniform prior
            for i in range(1,self.n_params):
                cube[i] = cube[i] * (self.X_up-self.X_low) + self.X_low
            # cube[1] = cube[1] * (self.X_up-self.X_low) + self.X_low
        

        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        DATASTDmean = mean(DATASTD)
        
        parameters = [str(i) for i in range(self.n_params)]
        n_params = self.n_params
        ndim = n_params




        # progress = pymultinest.ProgressPlotter(n_params = n_params); progress.start()
        # threading.Timer(60, show, ["chains/1-phys_live.points.pdf"]).start() # delayed opening
        pymultinest.run(multinest_loglike, multinest_uniform_prior, n_params,
                        importance_nested_sampling = self.params.nest_imp_sampling, resume = resume,
                        verbose = self.params.nest_verbose,sampling_efficiency = self.params.nest_samp_eff,
                        n_live_points = self.params.nest_nlive,max_iter= self.params.nest_max_iter,init_MPI=False)
        # progress.stop()

        #forcing slave processes to exit at this stage
        if MPI.COMM_WORLD.Get_rank() != 0:
            exit()
        
        
        #coallating results into arrays
        OUT = pymultinest.Analyzer(n_params=n_params)
        STATS = OUT.get_stats()
        Tout_mean, Tout_std, Xout_mean, Xout_std = self.collate_multinest_result(OUT)

        #saving arrays to object
        self.NEST        = True
        self.NEST_T_mean = Tout_mean
        self.NEST_T_std  = Tout_std
        self.NEST_X_mean = Xout_mean
        self.NEST_X_std  = Xout_std
        self.NEST_STATS  = STATS
        self.NEST_FITDATA= OUT


        
#     def multinest_uniform_prior(self,cube,ndim,nparams):
#         #prior distributions called by multinest
#         #implements a uniform prior
#
#         #converting temperatures from normalised grid to uniform prior
#         cube[0] = cube[0]* (self.T_up-self.T_low)+self.T_low
# #         print cube[0]
#
#         #converting mixing ratios from normalised grid to uniform prior
# #         for i in range(1,self.n_params):
# #             cube[i] = cube[i] * (self.X_up-self.X_low) + self.X_low
#
#         cube[1] = cube[1] * (self.X_up-self.X_low) + self.X_low
#
#         return cube

#     def multinest_loglike(self, cube, ndim,nparams):
#         #log-likelihood function called by multinest
#
#         PFIT = [cube[i] for i in range(self.n_params)]
#         PFIT = asarray(PFIT)
#
#         chi_t = self.chisq_trans(PFIT,DATA,DATASTD)
#         llterms =   -0.5* chi_t
#         return llterms
    