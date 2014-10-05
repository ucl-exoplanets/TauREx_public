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
import numpy as np
import os
import pylab as py
# import classes.emission as emi
from classes.transmission import *
from classes.profile import *
from library.library_plotting import *


class output(object):
    def __init__(self,params,data,fit=None):

        self.params = params
        self.data   = data
        self.fit    = fit
        
        #dumping fitted paramter list to file 
        with open('chains/parameters.dat','w') as parfile:
            parfile.write('Temperature \n')
            for mol in params.planet_molec:
                parfile.write(mol+' \n')
                
        
        #determining which fits were performed
        if fit is not None:
            self.MCMC = fit.MCMC
            self.NEST = fit.NEST
            self.DOWN = fit.DOWNHILL

        #calculating final absorption/emission spectrum (currently only absorption)
        self.profile = profile(params,data)
        if self.params.fit_transmission:
            self.trans = transmission(params,data,self.profile)
            if self.MCMC:
                self.transspec_mcmc = self.trans.cpath_integral(rho=self.profile.get_rho(T=fit.MCMC_T_mean),X=fit.MCMC_X_mean)
                self.transspec_mcmc = np.interp(self.data.wavegrid,self.data.specgrid,self.transspec_mcmc)
            if self.NEST:
                self.transspec_nest = self.trans.cpath_integral(rho=self.profile.get_rho(T=fit.NEST_T_mean),X=fit.NEST_X_mean)
                self.transspec_nest = np.interp(self.data.wavegrid,self.data.specgrid,self.transspec_nest)
                # print shape(fit.NEST_T_mean), shape(fit.NEST_X_mean)
            if self.DOWN:
                # print shape(fit.DOWNHILL_T_mean), shape(fit.DOWNHILL_X_mean)
                self.transspec_down = self.trans.cpath_integral(rho=self.profile.get_rho(T=fit.DOWNHILL_T_mean),X=fit.DOWNHILL_X_mean)
                self.transspec_down = np.interp(self.data.wavegrid,self.data.specgrid,self.transspec_down)
        elif self.params.fit_emission:
            print 'Emission not implemented yet'

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
                plot_mcmc_results(self.fit.MCMC_FITDATA,parameters=param_names,save2pdf=save2pdf, out_path=self.params.out_path)
            else:
                plot_mcmc_results(self.fit.MCMC_FITDATA,save2pdf=save2pdf, out_path=self.params.out_path)
        else:
            print 'MCMC was not run. Nothing to see here. Go away... '


    def plot_multinest(self,param_names=False,save2pdf=False):
        #plotting nested sampling distributions
        if self.NEST:
            if not param_names:
                parameters = range(self.fit.n_params)
            else:
                parameters = param_names
            plot_multinest_results(self.fit.NEST_FITDATA,parameters=parameters,save2pdf=save2pdf, out_path=self.params.out_path)
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
        py.ylabel('(Rp/Rs)^2') #dynamic for emission case later
        
        if save2pdf: fig.savefig(self.params.out_path+'spectrum+data.pdf')

    def plot_fit(self,save2pdf=False,linewidth=2.0):
        #plotting data and fits

        fig = py.figure()
        py.errorbar(self.data.spectrum[:,0],self.data.spectrum[:,1],self.data.spectrum[:,2],color=[0.7,0.7,0.7],linewidth=linewidth)
        py.plot(self.data.spectrum[:,0],self.data.spectrum[:,1],color=[0.3,0.3,0.3],linewidth=linewidth,label='DATA')
        if self.DOWN:
            py.plot(self.data.spectrum[:,0],transpose(self.transspec_down),c='b',label='DOWNHILL',linewidth=linewidth)
        if self.MCMC:
            py.plot(self.data.spectrum[:,0],transpose(self.transspec_mcmc),c='r',label='MCMC',linewidth=linewidth)
        if self.NEST:
            py.plot(self.data.spectrum[:,0],transpose(self.transspec_nest),c='g',label='NESTED',linewidth=linewidth)
        py.legend()
        py.title('Data and Model')
        py.xlabel('Wavelength ($\mu m$)')
        py.ylabel('(Rp/Rs)^2') #dynamic for emission case later
        
        if save2pdf: fig.savefig(self.params.out_path+'model_fit.pdf')
        

    def plot_manual(self,model,save2pdf=False,linewidth=2.0):
        fig = py.figure()
        py.plot(model[:,0],model[:,1],color=[0.0,0.0,0.0],linewidth=linewidth,label='Model')
        py.legend()
        py.title('Model')
        py.xlabel('Wavelength ($\mu m$)')
        py.ylabel('$(Rp/Rs)^2$') #dynamic for emission case later
        
        if save2pdf: fig.savefig(self.params.out_path+'spectrum.pdf')

    def save_model(self,ascii=True,modelout=None, modelsaveas= None):
        #dumps all model fits to file
        if not os.path.exists(self.params.out_path):
            os.makedirs(self.params.out_path)
        if self.params.fit_transmission:
            if self.MCMC and ascii:
                OUT = np.zeros((len(self.data.spectrum[:,0]),2))
                OUT[:,0] = self.data.spectrum[:,0]
                OUT[:,1] = np.transpose(self.transspec_mcmc)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+'transmission_spectrum_mcmc.dat',OUT)
            if self.NEST and ascii:
                OUT = np.zeros((len(self.data.spectrum[:,0]),2))
                OUT[:,0] = self.data.spectrum[:,0]
                OUT[:,1] = np.transpose(self.transspec_nest)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+'transmission_spectrum_nest.dat',OUT)
            if self.DOWN and ascii:
                OUT = np.zeros((len(self.data.spectrum[:,0]),2))
                OUT[:,0] = self.data.spectrum[:,0]
                OUT[:,1] = np.transpose(self.transspec_down)
                np.savetxt(self.params.out_path+self.params.out_file_prefix+'transmission_spectrum_down.dat',OUT)
        if self.params.fit_emission:
            print 'NEIN!'
        if modelout is not None:
            np.savetxt(self.params.out_path+modelsaveas,modelout)
            
            
