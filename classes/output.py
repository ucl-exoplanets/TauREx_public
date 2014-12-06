################################################
#class output
#
# Collates all modelling outputs and provides simple ascii output files
# also handles all plotting routines.
#
# Input: -parameter object
#        -data object
#        -fitting object
#
#
# Output: - plots and human readable files
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Jan 2014
#
################################################

#loading libraries
from base import base   
import numpy as np
import os
import pylab as py
import emission, transmission, tp_profile, library_plotting
from emission import *
from transmission import *
from tp_profile import *
from library_plotting import *


class output(base):
    def __init__(self, fitting, data=None, params=None):

        logging.info('Initialise object output')

        if params:
            self.params = params
        else:
            self.params = fitting.params #get params object from profile

        if data:
            self.data = data
        else:
            self.data = fitting.data    # get data object from profile

        self.fitting = fitting
        #self.profile  = tp_profile(params, data) #loading profile class (not needed, profile object is in fitting object)
        self.profile = fitting.profile
        
        #inheriting emission/transmission model from fitting object if loaded
        # 'fit' can either be the fitting instance or the emission/transmission instances
        
        if fitting.__ID__ == 'fitting':
            self.model        = fitting.model #directly inheriting model instance from fitting instance
            self.__MODEL_ID__ = fitting.__MODEL_ID__
        else:
            self.set_model(fit) #abstracting emission/transmission model 
        

        #dumping fitted paramter list to file 
        with open('chains/parameters.dat','w') as parfile:
            parfile.write('Temperature \n')
            for mol in self.params.planet_molec:
                parfile.write(mol+' \n')
                
        
        #determining which fits were performed
        if fitting.__ID__ == 'fitting':
            self.MCMC = fitting.MCMC
            self.NEST = fitting.NEST
            self.DOWN = fitting.DOWNHILL
        else:
            self.MCMC = False
            self.NEST = False
            self.DOWN = False

        #calculating final absorption/emission spectrum (currently only absorption)

        if self.MCMC:
            self.spec_mcmc = self.model(rho=self.profile.get_rho(T=fitting.MCMC_T_mean),X=fitting.MCMC_X_mean,temperature=fitting.MCMC_T_mean)
            self.spec_mcmc = np.interp(self.data.wavegrid,self.data.specgrid,self.spec_mcmc)
        if self.NEST:
            self.spec_nest = self.model(rho=self.profile.get_rho(T=fitting.NEST_T_mean),X=fitting.NEST_X_mean,temperature=fitting.NEST_T_mean)
            self.spec_nest = np.interp(self.data.wavegrid,self.data.specgrid,self.spec_nest)
        if self.DOWN:
            self.spec_down = self.model(rho=self.profile.get_rho(T=fitting.DOWNHILL_T_mean),X=fitting.DOWNHILL_X_mean,temperature=fitting.DOWNHILL_T_mean)
            self.spec_down = np.interp(self.data.wavegrid,self.data.specgrid,self.spec_down)

    #class methods
    


    def plot_all(self,save2pdf=False):
    #plots absolutely everything

        self.plot_spectrum(save2pdf=save2pdf)
        self.plot_fit(save2pdf=save2pdf)
        if self.MCMC:
            self.plot_mcmc(save2pdf=save2pdf)
        if self.NEST:
            self.plot_multinest(save2pdf=save2pdf)


    def plot_mcmc(self,param_names=False,save2pdf=False):
        #plotting MCMC sampled distributions
        if self.MCMC:
            if param_names:
                plot_mcmc_results(self.fitting.MCMC_FITDATA,parameters=param_names,save2pdf=save2pdf, out_path=self.params.out_path)
            else:
                plot_mcmc_results(self.fitting.MCMC_FITDATA,save2pdf=save2pdf, out_path=self.params.out_path)
        else:
            print 'MCMC was not run. Nothing to see here. Go away... '


    def plot_multinest(self,param_names=False,save2pdf=False):
        #plotting nested sampling distributions
        if self.NEST:
            if not param_names:
                parameters = range(self.fitting.n_params)
            else:
                parameters = param_names
            
            if self.params.nest_multimodes:
                plot_multinest_results(self.fitting.NEST_FITDATA,parameters=parameters,save2pdf=save2pdf, out_path=self.params.out_path)
            else:
                if self.params.verbose:
                    print 'WARNING: plotting routine for multimodes = False disabled. Please use plot_chains.py script.'

        else:
            print 'Multinest was not run. Nothing to see here. Go away... '


    def plot_spectrum(self,save2pdf=False,linewidth=2.0):
        #plotting data only

        fig = py.figure()
        py.errorbar(self.data.spectrum[:,0],self.data.spectrum[:,1],self.data.spectrum[:,2],color=[0.7,0.7,0.7],linewidth=linewidth)
        py.plot(self.data.spectrum[:,0],self.data.spectrum[:,1],color=[0.3,0.3,0.3],linewidth=linewidth,label='DATA')
        py.legend()
        py.title('Input data')
        py.xlabel('Wavelength')
        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$') 
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$') 
        
        if save2pdf: fig.savefig(self.params.out_path+'spectrum+data.pdf')

    def plot_fit(self,save2pdf=False,linewidth=2.0):
        #plotting data and fits

        fig = py.figure()
        py.errorbar(self.data.spectrum[:,0],self.data.spectrum[:,1],self.data.spectrum[:,2],color=[0.7,0.7,0.7],linewidth=linewidth)
        py.plot(self.data.spectrum[:,0],self.data.spectrum[:,1],color=[0.3,0.3,0.3],linewidth=linewidth,label='DATA')
        if self.DOWN:
            py.plot(self.data.spectrum[:,0],transpose(self.spec_down),c='b',label='DOWNHILL',linewidth=linewidth)
        if self.MCMC:
            py.plot(self.data.spectrum[:,0],transpose(self.spec_mcmc),c='r',label='MCMC',linewidth=linewidth)
        if self.NEST:
            py.plot(self.data.spectrum[:,0],transpose(self.spec_nest),c='g',label='NESTED',linewidth=linewidth)
        py.legend()
        py.title('Data and Model')
        py.xlabel('Wavelength ($\mu m$)')
        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$') 
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$') 
        
        if save2pdf: fig.savefig(self.params.out_path+'model_fit.pdf')
        

    def plot_manual(self,model,save2pdf=False,linewidth=2.0):
        fig = py.figure()
        py.plot(model[:,0],model[:,1],color=[0.0,0.0,0.0],linewidth=linewidth,label='Model')
        py.legend()
        py.title(self.__MODEL_ID__+' Model')
        py.xlabel('Wavelength ($\mu m$)')
        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$') 
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$') 
        
        if save2pdf: fig.savefig(self.params.out_path+'spectrum.pdf')

    def save_model(self,ascii=True,modelout=None, modelsaveas= None):
        #dumps all model fits to file
        if not os.path.exists(self.params.out_path):
            os.makedirs(self.params.out_path)
        
        OUT = np.zeros((len(self.data.spectrum[:,0]),2))
        OUT[:,0] = self.data.spectrum[:,0]
        
        if self.MCMC and ascii:
                OUT[:,1] = np.transpose(self.spec_mcmc)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+self.__MODEL_ID__+'_spectrum_mcmc.dat',OUT)
        if self.NEST and ascii:
                OUT[:,1] = np.transpose(self.spec_nest)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+self.__MODEL_ID__+'_spectrum_nest.dat',OUT)
        if self.DOWN and ascii:
                OUT[:,1] = np.transpose(self.spec_down)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+self.__MODEL_ID__+'_spectrum_down.dat',OUT)

        if modelout is not None:
            np.savetxt(self.params.out_path+modelsaveas,modelout)
            
            
