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
import emission, transmission, atmosphere, library_plotting
from emission import *
from transmission import *
from atmosphere import *
from library_plotting import *


class output(base):

    def __init__(self, fitting=None, forwardmodel=None, data=None, atmosphere=None, params=None):
         
        logging.info('Initialise object output')

        if params is not None:
            self.params = params
        else:
            if fitting:
                self.params = fitting.params #get params object from profile
            elif forwardmodel:
                self.params = forwardmodel.params

        if data is not None:
            self.data = data
        else:
            if fitting:
                self.data = fitting.data    # get data object from profile
            elif forwardmodel:
                self.data = forwardmodel.data

        if atmosphere is not None:
            self.atmosphere = atmosphere
        else:
            if fitting:
                self.atmosphere = fitting.atmosphere    # get data object from profile
            elif forwardmodel:
                self.atmosphere = forwardmodel.atmosphere

        if forwardmodel is not None:
            self.forwardmodel = forwardmodel
        else:
            self.forwardmodel = fitting.forwardmodel    # get data object from profile

        self.fitting = fitting

        self.__MODEL_ID__ = type(self.forwardmodel).__name__

        #dumping fitted paramter names to file and list

        self.parameters = []
        with open(os.path.join(self.params.out_path, 'parameters.dat'), 'w') as parfile:
            if not self.params.fit_fix_temp:
                parfile.write('Temperature \n')
                self.parameters.append('Temperature')
            if not self.params.fit_couple_mu:
                if not self.params.fit_fix_mu:
                    parfile.write('mu \n')
                    self.parameters.append('mu')
            if not self.params.fit_fix_radius:
                parfile.write('Radius \n')
                self.parameters.append('Radius')
            if not self.params.fit_fix_P0:
                parfile.write('P0 \n')
                self.parameters.append('P0')

            for mol in self.params.planet_molec:
                parfile.write('%s \n' % mol)
                self.parameters.append(mol)

            for mol in self.params.planet_inactive_gases:
                parfile.write('%s \n' % mol)
                self.parameters.append(mol)

        #determining which fits were performed
        if fitting is not None and fitting.__ID__ == 'fitting':
            self.MCMC = fitting.MCMC
            self.NEST = fitting.NEST
            self.DOWN = fitting.DOWNHILL
        else:
            self.MCMC = False
            self.NEST = False
            self.DOWN = False

        #calculating final absorption/emission spectrum (currently only absorption)
        # @todo Why reinterpolation to observed wavelengths?? We should probably bin to obs spectrum rather than interpolate?

        if self.MCMC:
            # @todo This is only for thread 0. Maybe output results of all chains, stored in fitting.MCMC_OUT[MPIRank]
            self.spec_mcmc = self.fitting.forwardmodel(rho=self.atmosphere.get_rho(T=self.fitting.MCMC_T_mean),
                                                X=fitting.MCMC_X_mean,
                                                temperature=self.fitting.MCMC_T_mean)
            self.spec_mcmc = np.interp(self.data.wavegrid, self.data.specgrid, self.spec_mcmc)

        if self.NEST:

            # separate solutions (modes) have different spectra. These are stored in self.spec_nest[i],

            self.spec_nest = []
            for fit_params in self.fitting.multinest_result:
                self.fitting.update_atmospheric_parameters(fit_params)
                self.spec_nest.append(self.fitting.forwardmodel.model())

            # interpolate
            self.spec_nest = [np.interp(self.data.wavegrid, self.data.specgrid, self.spec_nest[i]) \
                                        for i in range(self.fitting.NEST_modes)]

        if self.DOWN:
            self.spec_down = self.fitting.forwardmodel.model(rho=self.atmosphere.get_rho(T=self.fitting.DOWNHILL_T_mean),
                                                    X=self.fitting.DOWNHILL_X_mean,
                                                    temperature=self.fitting.DOWNHILL_T_mean)
            self.spec_down = np.interp(self.data.wavegrid,self.data.specgrid,self.spec_down)

    #class methods

    def plot_all(self, save2pdf=False):

        logging.info('Plotting absolutely everything')

        self.plot_spectrum(save2pdf=save2pdf)
        self.plot_fit(save2pdf=save2pdf)

        if self.MCMC:
            self.plot_mcmc(save2pdf=save2pdf, param_names=self.parameters)
        if self.NEST:
            self.plot_multinest(save2pdf=save2pdf, param_names=self.parameters)

    def plot_mcmc(self, param_names=False, save2pdf=False):

        logging.info('Plotting nested sampling distributions. Saving to %s' % self.params.out_path)

        if self.MCMC:

            if not param_names:
                parameters = range(self.fitting.n_params)
            else:
                parameters = param_names

            if param_names:
                plot_mcmc_results(self.fitting.dir_mcmc,
                                  parameters=parameters,
                                  save2pdf=save2pdf,
                                  out_path=self.params.out_path,
                                  plot_contour=self.params.out_plot_contour)
        else:
            logging.warning('MCMC was not run. Nothing to see here. Go away...')

    def plot_multinest(self, param_names=False, save2pdf=False):
        #plotting nested sampling distributions.
        # @todo this function is almost identical to plot_mcmc... Maybe merge?

        if self.NEST:

            if not param_names:
                parameters = range(self.fitting.n_params)
            else:
                parameters = param_names

            if self.params.nest_multimodes:

                logging.info('Plotting nested sampling distributions. Saving to %s' % self.params.out_path)

                plot_multinest_results(self.fitting.dir_multinest,
                                       parameters=self.fitting.fit_params_names,
                                       save2pdf=save2pdf,
                                       out_path=self.params.out_path,
                                       plot_contour=self.params.out_plot_contour)
            else:
                logging.warning('WARNING: plotting routine for multimodes = False disabled. '
                                'Please use plot_chains.py script.')

        else:
            logging.warning('MULTINEST was not run. Nothing to see here. Go away...')

    def plot_spectrum(self,save2pdf=False,linewidth=2.0):

        logging.info('Plotting observed spectrum')

        fig = py.figure()
        py.errorbar(self.data.spectrum[:,0],
                    self.data.spectrum[:,1],
                    self.data.spectrum[:,2],
                    color=[0.7,0.7,0.7],
                    linewidth=linewidth)

        # @todo why again?
        py.plot(self.data.spectrum[:,0],
                self.data.spectrum[:,1],
                color=[0.3,0.3,0.3],
                linewidth=linewidth,
                label='Data')

        py.legend()
        py.title('Input data')
        py.xlabel('Wavelength ($\mu m$)')

        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$')

        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$')

        if save2pdf:
            filename = os.path.join(self.params.out_path, 'spectrum_data.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)

    def plot_fit(self,save2pdf=False,linewidth=2.0):

        logging.info('Plotting observed and fitted spectra')

        fig = py.figure()

        # plot observed spectrum
        py.errorbar(self.data.spectrum[:,0],
                    self.data.spectrum[:,1],
                    self.data.spectrum[:,2],
                    color=[0.7,0.7,0.7],linewidth=linewidth)
        py.plot(self.data.spectrum[:,0],
                self.data.spectrum[:,1],
                color=[0.3,0.3,0.3],
                linewidth=linewidth,label='DATA')

        # plot models
        if self.DOWN:
            py.plot(self.data.spectrum[:,0],transpose(self.spec_down),c='b',label='DOWNHILL',linewidth=linewidth)
        if self.MCMC:
            py.plot(self.data.spectrum[:,0],transpose(self.spec_mcmc),c='r',label='MCMC',linewidth=linewidth)
        if self.NEST:
            for i in range(self.fitting.NEST_modes):
                py.plot(self.data.spectrum[:,0],transpose(self.spec_nest[i]),
                        label='NESTED %i' % i, linewidth=linewidth)

        py.legend()
        py.title('Data and Model')
        py.xlabel('Wavelength ($\mu m$)')
        py.xscale('log')
        py.xlim((np.min(self.data.spectrum[:,0]), np.max(self.data.spectrum[:,0])))
        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$')
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$')

        if save2pdf:
            filename = os.path.join(self.params.out_path, 'model_fit.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)


    def plot_manual(self,model,save2pdf=False,linewidth=2.0):

        fig = py.figure()
        py.plot(model[:,0],model[:,1],color=[0.0,0.0,0.0],linewidth=linewidth,label='Model')
        py.legend()
        py.xscale('log')
        py.title(self.__MODEL_ID__+' Model')
        py.xlabel('Wavelength ($\mu m$)')

        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$')
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$')

        if save2pdf:
            filename = os.path.join(self.params.out_path, 'spectrum.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)


    def save_model(self, ascii=True, modelout=None, modelsaveas=None):

        logging.info('Dumps all model fits to file')

        if modelout == None:
            out = np.zeros((len(self.data.spectrum[:,0]), 2))
            out[:,0] = self.data.spectrum[:,0]
        else:
            out = modelout

        basename = os.path.join(self.params.out_path, self.params.out_file_prefix + self.__MODEL_ID__)

        if self.MCMC and ascii:
                out[:,1] = np.transpose(self.spec_mcmc)
                filename = str(basename)+'_spectrum_mcmc.dat'
                logging.info('Saving MCMC spectrum to %s' % filename)
                np.savetxt(filename, out)

        if self.NEST and ascii:
            logging.info('MultiNest detected %i different modes. '
                         'Saving one model spectrum for each solution' % self.fitting.NEST_modes)
            for i in range(self.fitting.NEST_modes):
                out[:,1] = np.transpose(self.spec_nest[i])
                filename = '%s_spectrum_nest_%i.dat' % (basename, i)
                logging.info('Saving Nested Sampling spectrum for solution %i to %s' % (i, filename))
                np.savetxt(filename, out)

        if self.DOWN and ascii:
                out[:,1] = np.transpose(self.spec_down)
                filename = str(basename)+'_spectrum_down.dat'
                logging.info('Saving MLE spectrum to %s' % filename)
                np.savetxt(filename, out)

        if modelout is not None: # ???
            np.savetxt(self.params.out_path+modelsaveas,modelout)



