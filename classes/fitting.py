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
import numpy, pylab
from numpy import *
from pylab import *
from StringIO import StringIO
from scipy.interpolate import interp1d
from scipy.optimize import fmin
import pymc
import warnings


class fitting(object):
    def __init__(self,params,data,profile,transmod):
        
        self.params      = params
        self.dataob      = data
        self.transmod    = transmod
        self.profileob   = profile
        self.observation = data.spectrum
        self.nlayers = params.tp_atm_levels
        self.ngas    = params.tp_num_gas

              
        
        #checking if number of free parameters for X = number of gas species
#         if len(self.params.fit_param_free_X) != (self.params.tp_atm_levels*self.params.tp_num_gas):
#             warnings.warn('Number of free X parameters does not equal number of gas species and atmospheric levels')
#             self.params.fit_param_free_X = self.params.tp_atm_levels*self.params.tp_num_gas
        
        #setting up initial parameter array
        self.PINIT = zeros((self.nlayers,self.ngas+1))
        self.PINIT[:,0] = self.profileob.T
        self.PINIT[:,1:] = self.profileob.X
        
        self.Pshape = shape(self.PINIT)
        
        self.PINIT2 = self.PINIT.flatten(order='F')
        self.PINIT3 = self.PINIT2[(self.Pshape[0]-1):] #restricting temperature to one free parameter only

       
#         figure(1)
#         plot(self.observation[:,1])
#         plot(self.transmod.cpath_integral(rho=self.profileob.get_rho(T=100)),'r')
#         show()
#         
#         exit()
        
#         rho = self.profileob.get_rho(T=self.PINIT[:,0]) 
#         print rho
#         exit()
        
        
    def chisq_trans(self,PFIT,DATA,DATASTD):       
        
        rho = self.profileob.get_rho(T=PFIT[0]) 
        X   = PFIT[1:].reshape(self.Pshape[0],self.Pshape[1]-1)
        
        MODEL = self.transmod.cpath_integral(rho=rho,X=X)
        
        res = (DATA-MODEL) / DATASTD
        return sum(res**2)
        
        
        
        
    def downhill_fit(self):
    # fits a Mandel & Agol lightcurve using simplex downhill minimisation

        PINIT = self.PINIT3  # initial temperatures and abundances
        
        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        
        PFIT, err, out3, out4, out5 = fmin(self.chisq_trans, PINIT, args=(DATA,DATASTD), disp=1, full_output=1)

        return PFIT
    
    
    def mcmc_fit(self):
    # adaptive Markov Chain Monte Carlo
    
        PINIT = self.PINIT3
        DATA = self.observation[:,1] #observed data
        DATASTD = self.observation[:,2] #data error 
        
        # setting prior distributions
        priors = empty(size(PINIT),dtype=object)
        priors[0] = pymc.Uniform('temp_prior',self.params.planet_temp-300.0,self.params.planet_temp+300.0) #uniform temperature prior
        for i in range(1,size(PINIT)):
            priors[i] = pymc.Uniform('mixing_prior_%i' % i, 0.0,1e-4)  # uniform mixing ratio prior 
        

        #setting up data error prior if specified
        if self.params.mcmc_update_std:
            std_dev = pymc.Uniform('std_dev',0.0,max(DATASTD),value=mean(DATASTD),N=len(DATASTD)) #uniform prior on data standard deviation
        else:
            std_dev = pymc.Uniform('std_dev',0.0,max(DATASTD),value=DATASTD,observed=True,N=len(DATASTD))
        
        
        @pymc.deterministic(plot=False)
        def precision(std_dev=std_dev):
                return std_dev
    
        
        # log-likelihood function
        @pymc.stochastic(observed=True, plot=False)
        def mcmc_loglikelihood(value=PINIT, PFIT=priors, DATASTD=precision, DATA=DATA):
             
            chi_t = self.chisq_trans(PFIT,DATASTD,DATA)
            # llterms = len(DATA) * log(1.0 / sqrt(2.0 * pi)) - len(DATA) * log(DATASTD) - 0.5 * chi_t
            llterms =  - 0.5* chi_t
            return llterms

        # executing MCMC sampling
        R = pymc.MCMC((priors, mcmc_loglikelihood), verbose=2)  # build the model
        R.sample(iter=self.params.mcmc_iter, burn=self.params.mcmc_burn,
                 thin=self.params.mcmc_thin)              # populate and run it

    
        return R