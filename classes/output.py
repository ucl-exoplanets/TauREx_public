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
import pylab as py
import os, glob, pickle, shutil
import emission, transmission, atmosphere, library_plotting, library_emission
from emission import *
from transmission import *
from atmosphere import *
from library_plotting import *
from library_emission import *
from scipy.optimize import curve_fit

#conversion constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg


class output(base):


    def __init__(self, fitting=None, forwardmodel=None, data=None, atmosphere=None, params=None, out_path=None):

        logging.info('Initialise object output')
        
        if params is not None:
            self.params = params
        else:
            if fitting:
                self.params = fitting.params #get params object from profile
            elif forwardmodel:
                self.params = forwardmodel.params
                
        if out_path is None:
            self.out_path = self.params.out_path
        else:
            self.out_path = out_path

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

        if fitting is not None:
            self.fitting = fitting

        if forwardmodel is not None:
            self.forwardmodel = forwardmodel
        else:
            self.forwardmodel = fitting.forwardmodel    # get data object from profile

        self.__MODEL_ID__ = type(self.forwardmodel).__name__


        if fitting is not None:
            self.NEST = self.fitting.NEST
            self.MCMC = self.fitting.MCMC
            self.DOWN = self.fitting.DOWN
        else:
            self.NEST = False
            self.MCMC = False
            self.DOWN = False

        if self.DOWN:
            self.store_downhill_solution()
        if self.MCMC:
            self.store_mcmc_solutions()
        if self.NEST:
            self.store_nest_solutions()

        if fitting is not None:
            self.save_fit_output_files()

    #class methods

    def store_downhill_solution(self):

        logging.info('Store the DOWNHILL results')

        self.params_names = self.fitting.fit_params_names

        dict = {}
        for idx, param_name in enumerate(self.params_names):
            dict[param_name] = {'value': self.fitting.DOWN_fit_output[idx] }
        DOWN_out = [{'fit_params': dict}]

        self.DOWN_out = self.add_spectra_from_solutions(DOWN_out)
        self.DOWN_params_values = self.fitting.DOWN_fit_output
        self.DOWN_TP_params_values = self.DOWN_params_values[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fitting.fit_TP_nparams]



    def store_mcmc_solutions(self):

        logging.info('Store the MCMC results')

        # read individual traces from PyMCMC output and combine them into a single array, called tracedata
        for thread in glob.glob(os.path.join(self.fitting.dir_mcmc, '*')):
            chainlist = glob.glob(os.path.join(thread, 'Chain_0', 'PFIT*'))
            tracedata = [np.loadtxt(trace) for trace in chainlist]
        tracedata = np.transpose(tracedata)

        if self.params.fit_clr_trans:
            # the mixing ratios are sampled in the log-ratio space. Need to convert them back to simplex
            self.MCMC_out, self.MCMC_tracedata_clr_inv, self.MCMC_labels =  self.analyse_traces(tracedata, clr_inv=True)
            self.MCMC_tracedata = tracedata
        else:
            self.MCMC_out, self.MCMC_tracedata, self.MCMC_labels = self.analyse_traces(tracedata)

        self.MCMC_params_values = [self.MCMC_out[0]['fit_params'][param]['value'] for param in self.params_names]
        self.MCMC_params_std = [self.MCMC_out[0]['fit_params'][param]['std'] for param in self.params_names]
        self.MCMC_X_params_values = self.MCMC_params_values[:self.fitting.fit_X_nparams]
        self.MCMC_X_params_std = self.MCMC_params_std[:self.fitting.fit_X_nparams]
        self.MCMC_TP_params_values = self.MCMC_params_values[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fitting.fit_TP_nparams]
        self.MCMC_TP_params_std = self.MCMC_params_std[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fitting.fit_TP_nparams]

    def store_nest_solutions(self):

        logging.info('Store the MULTINEST results')

        # get output traces from 1-.txt. Exclude first two columns.
        NEST_data = np.loadtxt(os.path.join(self.fitting.dir_multinest, '1-.txt'))

        self.NEST_tracedata = NEST_data[:,2:] # exclude first two columns
        self.NEST_likelihood = NEST_data[:,:2] # get first two columns

        self.NEST_out, self.NEST_tracedata =  self.get_multinest_solutions()

        # storing fit_params output and standard deviation
        self.NEST_params_values    = {}
        self.NEST_params_std       = {}
        self.NEST_TP_params_std    = {}
        self.NEST_X_params_values  = {}
        self.NEST_X_params_std     = {}
        self.NEST_TP_params_values = {}
        self.NEST_TP_params_std    = {}
        for idx, solution in enumerate(self.NEST_out):
            self.NEST_params_values[idx] = [self.NEST_out[idx]['fit_params'][param]['value'] for param in self.params_names]
            self.NEST_params_std[idx]    = [self.NEST_out[idx]['fit_params'][param]['std'] for param in self.params_names]
            self.NEST_X_params_values[idx]  = self.NEST_params_values[idx][:self.fitting.fit_X_nparams]
            self.NEST_X_params_std[idx]     = self.NEST_params_std[idx][:self.fitting.fit_X_nparams]
            self.NEST_TP_params_values[idx] = self.NEST_params_values[idx][self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fitting.fit_TP_nparams]
            self.NEST_TP_params_std[idx]    = self.NEST_params_std[idx][self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fitting.fit_TP_nparams]

    def get_multinest_solutions(self):

        # get solutions from multinest
        self.params_names = self.fitting.fit_params_names

        data = np.loadtxt(os.path.join(self.fitting.dir_multinest, '1-.txt'))

        NEST_analyzer = pymultinest.Analyzer(n_params=len(self.fitting.fit_params),
                                             outputfiles_basename=os.path.join(self.fitting.dir_multinest, '1-'))
        NEST_stats = NEST_analyzer.get_stats()
        NEST_out = []

        modes = []
        modes_weights = []
        chains = []
        chains_weights = []
        if self.params.nest_multimodes:
            # get parameter values and sample probability (=weight) for each mode
            with open(os.path.join(self.fitting.dir_multinest, '1-post_separate.dat')) as f:
                lines = f.readlines()
                for idx, line in enumerate(lines):
                    if idx > 2: # skip the first two lines
                        if lines[idx-1] == '\n' and lines[idx-2] == '\n':
                            modes.append(chains)
                            modes_weights.append(chains_weights)
                            chains = []
                            chains_weights = []
                    chain = [float(x) for x in line.split()[2:]]
                    if len(chain) > 0:
                        chains.append(chain)
                        chains_weights.append(float(line.split()[0]))
                modes.append(chains)
                modes_weights.append(chains_weights)
            modes_array = []
            for mode in modes:
                mode_array = np.zeros((len(mode), len(mode[0])))
                for idx, line in enumerate(mode):
                    mode_array[idx, :] = line
                modes_array.append(mode_array)
        else:
            # not running in multimode. Get chains directly from file 1-.txt
            modes_array = [data[:,2:]]
            chains_weights = [data[:,0]]

        NEST_tracedata = data[:,2:]

        for nmode in range(len(modes)):

            dict = {'weights': modes_weights[nmode],
                    'fit_params': {}}
            for idx, param_name in enumerate(self.params_names):
                dict['fit_params'][param_name] = {
                    'value': NEST_stats['modes'][nmode]['maximum a posterior'][idx],
                    'std': NEST_stats['modes'][nmode]['sigma'][idx],
                    'trace': modes_array[nmode][:,idx],
                }

            if self.params.fit_clr_trans == True:

                nallgases = self.fitting.forwardmodel.atmosphere.nallgases

                # if centered-log-ratio transformation is True, transform the traces of the log ratios back to abundances
                mixing_ratios_clr = modes_array[nmode][:, :nallgases-1] # get traces of log-ratios
                mixing_ratios_clr = np.c_[mixing_ratios_clr, -sum(mixing_ratios_clr, axis=1)] #add last clr (= -sum(clr) )
                mixing_ratios_clr_inv, coupled_mu_trace = self.inverse_clr_transform(mixing_ratios_clr)

                self.NEST_data_clr_inv = np.c_[mixing_ratios_clr_inv, coupled_mu_trace]

                # get means
                clr = np.asarray([NEST_stats['modes'][nmode]['maximum a posterior'][i] for i in range(nallgases-1)])
                clr = np.append(clr, -np.sum(clr))

                # add last CLR to dictionary. This is not fitted, as it is equal to minus the sum of all the other CLRs
                # Note that that the sum of all CLR is equal to zero
                dict['fit_params']['%s_CLR' % self.forwardmodel.atmosphere.inactive_gases[-1]] = {
                    'value': clr[nallgases-1],
                    'std': weighted_avg_and_std(mixing_ratios_clr[:, nallgases-1], modes_weights[nmode])[1], # weighted std
                    'trace': mixing_ratios_clr[:, nallgases-1],
                }

                mixing_means = np.exp(clr) # add log-ratio (= -sum(other log ratios)
                mixing_means /= sum(mixing_means) # closure operation

                # get coupled std of means
                absorbing_gases_X = mixing_means[:len(self.fitting.forwardmodel.atmosphere.absorbing_gases)]
                inactive_gases_X = mixing_means[len(self.fitting.forwardmodel.atmosphere.inactive_gases)+1:]

                coupled_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu(absorbing_gases_X, inactive_gases_X)

                # convert mixing means to log space # todo consider new parameter X_log
                mixing_means = np.log10(mixing_means)

                # set new list of params names, including the clr inv mixing ratios
                self.clrinv_params_names = []
                for idx, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases +
                        self.fitting.forwardmodel.atmosphere.inactive_gases):
                    self.clrinv_params_names.append(gasname)
                self.clrinv_params_names.append('coupled_mu')
                for i in range(self.fitting.forwardmodel.atmosphere.nallgases-1, len(self.fitting.fit_params_names)):
                     self.clrinv_params_names.append(self.fitting.fit_params_names[i])

                for idx in range(self.fitting.forwardmodel.atmosphere.nallgases+1): # all gases + coupled_mu
                    if self.clrinv_params_names[idx] == 'coupled_mu':
                        trace = coupled_mu_trace
                        value = coupled_mu
                    else:
                        trace = mixing_ratios_clr_inv[:,idx]
                        value = mixing_means[idx]
                    std = weighted_avg_and_std(trace, modes_weights[nmode])[1]
                    dict['fit_params'][self.clrinv_params_names[idx]] = {
                        'value': value,
                        'std': std,
                        'trace': trace,
                    }

            NEST_out.append(dict)

        NEST_out = self.add_spectra_from_solutions(NEST_out)
        return NEST_out, NEST_tracedata

    def inverse_clr_transform(self, clr_tracedata):

        # Convert the traces of the log-ratios back to the gas mixing ratios in log space. This is done by calculating the inverse
        # CLR in the tracedata array, line by line. The new traces are then returned.

        logging.info('Inverse log-ratio transformation of the traces')

        ntraces = len(clr_tracedata[:,0])
        ncol = len(clr_tracedata[0,:])

        mixing_ratios_tracedata = exp(clr_tracedata)
        mixing_ratios_tracedata /= (np.asarray([np.sum(mixing_ratios_tracedata, axis=1),]*ncol).transpose()) # closure operation

        # calculate couple mean molecular weight trace as a byproduct
        coupled_mu_trace = np.zeros(ntraces)
        i = 0
        for sample in mixing_ratios_tracedata:
            # get mean molecular weight from gas mixing ratios
            absorbing_gases_X = sample[:len(self.fitting.forwardmodel.atmosphere.absorbing_gases)]
            inactive_gases_X = sample[len(self.fitting.forwardmodel.atmosphere.inactive_gases)+1:]
            coupled_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu(absorbing_gases_X, inactive_gases_X)
            coupled_mu_trace[i] = coupled_mu
            i += 1

        return np.log10(mixing_ratios_tracedata), coupled_mu_trace

    def add_spectra_from_solutions(self, fitting_out):

        # loop through solutions
        for idx, solution in enumerate(fitting_out):

            fit_params = [solution['fit_params'][param]['value'] for param in self.params_names]
            self.fitting.update_atmospheric_parameters(fit_params)
            model = self.fitting.forwardmodel.model()
            model_binned = [model[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid+1)]

            # model spectrum binned to observed spectrum
            out = np.zeros((len(self.data.spectrum[:,0]), 2))
            out[:,0] = self.data.spectrum[:,0]
            out[:,1] = model_binned
            fitting_out[idx]['spectrum'] = out

            # high res spectrum
            out = np.zeros((len(self.data.specgrid), 2))
            out[:,0] = self.data.specgrid
            out[:,1] = model
            fitting_out[idx]['highres_spectrum'] = out
            
            #individual molecules plotted @todo may need some cleaning 
            fit_params2 = fit_params           
            fitting_out[idx]['components'] ={}
            for idx2, param in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
                fit_params2 = np.copy(fit_params)
                if self.params.fit_X_log:
                    fit_params2[:self.fitting.fit_X_nparams] = [-20.0]*self.fitting.fit_X_nparams 
                else:
                    fit_params2[:self.fitting.fit_X_nparams] = [0.0]*self.fitting.fit_X_nparams
                fit_params2[idx2] = solution['fit_params'][param]['value']
#                 fit_params2[idx2] = 0.0
           
                
                self.fitting.update_atmospheric_parameters(fit_params2)
                model2 = self.fitting.forwardmodel.model()
                out2 = np.zeros((len(self.data.specgrid), 2))
                out2[:,0] = self.data.specgrid
                out2[:,1] = model2
                out3 = np.zeros((len(self.data.spectrum[:,0]), 2))
                out3[:,0] = self.data.spectrum[:,0]
                out3[:,1] = [model2[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid+1)]
                
                fitting_out[idx]['components'][param] ={}
                fitting_out[idx]['components'][param]['highres_spectrum'] = out2
                fitting_out[idx]['components'][param]['spectrum'] = out3

        return fitting_out

    def save_fit_output_files(self):

        # store parameter names of traces
        f = open(os.path.join(self.out_path, 'parameters.txt'),'w')
        for param in self.fitting.fit_params_names:
            f.write('%s\n' % param)
        f.close()

        # if clr transform, store parameter names of transformed traces
        if self.params.fit_clr_trans:

            try:
                self.params_names
                # If clr inverse, write out parameter names with mixing ratio naems
                f = open(os.path.join(self.out_path, 'parameters_clr_inv.txt'),'w')
                for param in self.params_names:
                    f.write('%s\n' % param)
                f.close()
            except: pass

        if self.DOWN:
            pickle.dump(self.DOWN_out, open(os.path.join(self.out_path, 'DOWN_out.db'), 'wb'))

        if self.MCMC:
            np.savetxt(os.path.join(self.out_path, 'MCMC_tracedata.txt'), self.MCMC_tracedata)
            if self.params.fit_clr_trans:
                np.savetxt(os.path.join(self.out_path, 'MCMC_clr_inv_tracedata.txt'), self.MCMC_data_clr_inv)
            pickle.dump(self.MCMC_out, open(os.path.join(self.out_path, 'MCMC_out.db'), 'wb'))
            self.save_fit_out_to_file(self.MCMC_out, type='MCMC')

        if self.NEST:
            np.savetxt(os.path.join(self.out_path, 'NEST_tracedata.txt'), self.NEST_tracedata)
            np.savetxt(os.path.join(self.out_path, 'NEST_likelihood.txt'), self.NEST_likelihood)
            if self.params.fit_clr_trans:
                np.savetxt(os.path.join(self.out_path, 'NEST_clr_inv_tracedata.txt'), self.NEST_data_clr_inv)
            pickle.dump(self.NEST_out, open(os.path.join(self.out_path, 'NEST_out.db'), 'wb'))
            self.save_fit_out_to_file(self.NEST_out, type='NEST')

        # save param file and observation to files
        try:
            shutil.copy(self.params.parfile, self.out_path)
            shutil.copy(self.params.in_spectrum_file, self.out_path)
        except:
            pass

    def save_fit_out_to_file(self, fit_out, type=''):

        ''' Save value and standard deviation of each parameter for each solution to txt file '''

        f = open(os.path.join(self.out_path, '%s_out.txt' % type),'w')
        f.write('Note: mixing ratios are expressed as fractions in linear space, mean molecular weight '
                'in AMU, pressure in Pascal\n\n')

        for idx, solution in enumerate(fit_out):
            f.write('%s solution %i\n' % (type, idx))

            for param in solution['fit_params'].keys():

                if param in self.params.all_absorbing_gases or param in self.params.all_inactive_gases:
                    if self.params.fit_X_log:
                        f.write('log(%s)	%f	%f\n' %  (param, solution['fit_params'][param]['value'],
                                                         solution['fit_params'][param]['std']))
                        f.write('%s	%f	%f\n' %  (param, np.power(10, solution['fit_params'][param]['value']),
                                                  np.power(10, solution['fit_params'][param]['value'])*np.log(10)*solution['fit_params'][param]['std'] ))
                    else:
                        f.write('%s	%f	%f\n' %  (param, solution['fit_params'][param]['value'],
                                                         solution['fit_params'][param]['std']))
                else:
                    if param == 'P0':
                        f.write('log(%s)	%f	%f\n' %  (param, solution['fit_params'][param]['value'],
                                                         solution['fit_params'][param]['std']))
                        f.write('%s	%f	%f\n' %  (param, np.power(10, solution['fit_params'][param]['value']),
                                                  np.power(10, solution['fit_params'][param]['value'])*np.log(10)*solution['fit_params'][param]['std'] ))
                    elif param == 'mu' or param == 'coupled_mu':
                        f.write('%s (AMU)	%f	%f\n' %  (param, solution['fit_params'][param]['value'],
                                                          solution['fit_params'][param]['std']))
                    else:
                        f.write('%s	%.3e	%.3e\n' %  (param, solution['fit_params'][param]['value'],
                                solution['fit_params'][param]['std']))

            f.write('\n')

    def save_fig(self,fig,FNAME):
            filename = os.path.join(self.out_path,FNAME)
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)

    def plot_all(self, save2pdf=False,params_names=None):

        logging.info('Plotting absolutely everything')

#         self.plot_spectrum(save2pdf=save2pdf)
        self.plot_fit(save2pdf=save2pdf,resolution='low')
        self.plot_fit(save2pdf=save2pdf,resolution='high')
        self.plot_absorbers(save2pdf=save2pdf,resolution='high')
        self.plot_distributions(save2pdf=save2pdf, params_names=params_names)

    def plot_spectrum(self,save2pdf=False,linewidth=2.0):

        logging.info('Plotting observed spectrum')

        fig = py.figure()
        py.errorbar(self.data.spectrum[:,0],
                    self.data.spectrum[:,1],
                    self.data.spectrum[:,2],
                    color=[0.7,0.7,0.7],
                    linewidth=linewidth)
        py.plot(self.data.spectrum[:,0],
                self.data.spectrum[:,1],
                color=[0.3,0.3,0.3],
                linewidth=linewidth,
                label='Data')
        py.xscale('log')


        py.legend()
        py.title('Input data')
        py.xlabel('Wavelength ($\mu m$)')
        py.xscale('log')
        py.xlim((np.min(self.data.spectrum[:,0]), np.max(self.data.spectrum[:,0])))

        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$')
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$')

        if save2pdf:
            self.save_fig(fig, 'spectrum_data.pdf')
        
            

    def plot_fit(self,save2pdf=False,resolution='low',linewidth=1.0):

        logging.info('Plotting observed and fitted spectra')

        if resolution is 'low':
            res = 'spectrum'
        elif resolution is 'high':
            res = 'highres_spectrum'

        fig = py.figure()

        # plot models       
        if self.fitting.DOWN:
            fig = self.__plot_fit__(self.DOWN_out[0][res], 'DOWNHILL', 
                                    fig=fig, linewidth=linewidth)
        if self.fitting.MCMC:
            for idx, solution in enumerate(self.MCMC_out):
                fig = self.__plot_fit__(solution[res], 'MCMC %i' % idx,
                                        fig=fig, linewidth=linewidth)
        if self.fitting.NEST:
            for idx, solution in enumerate(self.NEST_out):
                fig = self.__plot_fit__(solution[res], 'NESTED %i' % idx,
                                        fig=fig, linewidth=linewidth)

        # plot observed spectrum
        py.errorbar(self.data.spectrum[:,0],
                        self.data.spectrum[:,1],
                        self.data.spectrum[:,2],
                        color=[0.4,0.4,0.4], mec=[0.5,0.5,0.5], fmt='', marker='o', linewidth=2, linestyle='None')
        # py.plot(self.data.spectrum[:,0], # todo why again? - it looks prettier! Ah ok :)
        #             self.data.spectrum[:,1],
        #             color=[0.3,0.3,0.3],
        #             linewidth=0,label='DATA')

        if save2pdf:
            self.save_fig(fig,'model_fit_{0}_res.pdf'.format(resolution))
           


    def plot_absorbers(self,save2pdf=False,params_names=None,resolution='low',linewidth=2.0):
        #routine plotting individual components
        logging.info('Plotting individual absorbers.')
        if params_names is None:
            try:
                self.params_names
                params_names = self.params_names
            except:
                return
        
        if resolution is 'low':
            res = 'spectrum'
        elif resolution is 'high':
            res = 'highres_spectrum'
        
        
        # plot models
        if self.fitting.DOWN:
            fig = py.figure()
#             plot_observed = True
            fig = self.__plot_fit__(self.DOWN_out[0][res], 'MODEL', fig=fig, plot_observed = True,linewidth=linewidth)
            for idx2, param in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
#                 if idx2 == 1: plot_observed = False
                fig = self.__plot_fit__(self.DOWN_out[0]['components'][param][res], param, 
                                        fig=fig,plot_observed=False,linewidth=linewidth)
            if save2pdf:
                self.save_fig(fig, 'downhill_components_{}_res.pdf'.format(resolution))
            
        if self.fitting.MCMC:
            for idx, solution in enumerate(self.MCMC_out):
                fig = py.figure()
#                 plot_observed=True
                fig = self.__plot_fit__(solution[res], 'MODEL',fig=fig, plot_observed = True,linewidth=linewidth)
                for idx2, param in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
#                     if idx2 == 1: plot_observed = False
                    fig = self.__plot_fit__(solution['components'][param][res],param,
                                        fig=fig,plot_observed=False,linewidth=linewidth)
                    
                if save2pdf:
                    self.save_fig(fig, 'mcmc_components_{0}_{1}_res.pdf'.format(idx,resolution))
                
        if self.fitting.NEST:
            for idx, solution in enumerate(self.NEST_out):
                fig = py.figure()
#                 plot_observed=True
                fig = self.__plot_fit__(solution[res], 'MODEL',fig=fig, plot_observed = True,linewidth=linewidth)
                for idx2, param in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
#                     if idx2 == 1: plot_observed = False
                    fig = self.__plot_fit__(solution['components'][param][res],param,
                                        fig=fig,plot_observed=False,linewidth=linewidth)
                    
                if save2pdf:
                    self.save_fig(fig, 'nested_components_{0}_{1}_res.pdf'.format(idx,resolution))
                
    
    def __plot_fit__(self,MODEL,LABEL,fig=None,plot_observed=False, linewidth=2.0):
        
        if fig is None:
            fig = py.figure()

        if plot_observed:
            # plot observed spectrum
            py.errorbar(self.data.spectrum[:,0],
                        self.data.spectrum[:,1],
                        self.data.spectrum[:,2],
                        color=[0.7,0.7,0.7],linewidth=linewidth,label='DATA')
#             py.plot(self.data.spectrum[:,0], # todo why again? - it looks prettier!
#                     self.data.spectrum[:,1],
#                     color=[0.3,0.3,0.3],
#                     linewidth=linewidth,label='DATA')

        # plot models
        py.plot(MODEL[:,0], MODEL[:,1],label=LABEL,linewidth=linewidth, zorder=-32)
        
        py.legend(loc=0)
        py.title('Data and Model')
        py.xlabel('Wavelength ($\mu m$)')
        py.xscale('log')
        py.xlim((np.min(self.data.spectrum[:,0]), np.max(self.data.spectrum[:,0])))

        if self.__MODEL_ID__ == 'transmission':
            py.ylabel('$(Rp/Rs)^2$')
        elif self.__MODEL_ID__ == 'emission':
            py.ylabel('$F_p/F_s$')

        return fig
          

    def plot_distributions(self, save2pdf=False, params_names=None):

        logging.info('Plotting sampling distributions. Saving to %s' % self.out_path)

        if params_names is None:
            try:
                self.params_names
                params_names = self.params_names
            except:
                return

        if self.fitting.MCMC:
            plot_posteriors(self.MCMC_out,params_names=params_names,
                            save2pdf=save2pdf,out_path=self.out_path,plot_name = 'MCMC',
                            plot_contour=self.params.out_plot_contour,color=self.params.out_plot_colour)
        if self.fitting.NEST:
            plot_posteriors(self.NEST_out, params_names=params_names,save2pdf=save2pdf,out_path=self.out_path,
                            plot_name='NEST',plot_contour=self.params.out_plot_contour, color=self.params.out_plot_colour)

            if self.params.fit_clr_trans == True:
                params_names = self.clrinv_params_names
                plot_posteriors(self.NEST_out,params_names=params_names,save2pdf=save2pdf,out_path=self.out_path,
                                plot_name = 'NEST_clrinv',plot_contour=self.params.out_plot_contour, color=self.params.out_plot_colour)

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
            self.save_fig(fig, 'spectrum.pdf')


    def save_ascii_spectra(self):

        logging.info('Dumps all model fits and observed spectra to file')

        if self.MCMC:
            for idx, solution in enumerate(self.NEST_out):
                filename1 = os.path.join(self.out_path, 'MCMC_%s_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                filename2 = os.path.join(self.out_path, 'MCMC_%s_highres_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                logging.info('Saving MCMC spectrum for solution %i to %s and %s' % (idx, filename1, filename2))
                np.savetxt(filename1, solution['spectrum'])
                np.savetxt(filename2, solution['highres_spectrum'])

        if self.NEST:
            logging.info('MultiNest detected %i different modes. '
                         'Saving one model spectrum for each solution' % len(self.NEST_out))
            for idx, solution in enumerate(self.NEST_out):
                filename1 = os.path.join(self.out_path, 'NEST_%s_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                filename2 = os.path.join(self.out_path, 'NEST_%s_highres_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                logging.info('Saving Nested Sampling spectrum for solution %i to %s and %s' % (idx, filename1, filename2))
                np.savetxt(filename1, solution['spectrum'])
                np.savetxt(filename2, solution['highres_spectrum'])
        if self.DOWN:
                filename1 = os.path.join(self.out_path, 'DOWN_%s_spectrum.dat' % (self.__MODEL_ID__))
                filename2 = os.path.join(self.out_path, 'DOWN_%s_highres_spectrum.dat' % (self.__MODEL_ID__))
                logging.info('Saving DOWHNILL spectrum to %s and %s' % (filename1, filename2))
                np.savetxt(filename1, self.DOWN_out[0]['spectrum'])
                np.savetxt(filename2, self.DOWN_out[0]['highres_spectrum'])

        filename = os.path.join(self.out_path, 'observed_%s_spectrum.dat' % (self.__MODEL_ID__))
        np.savetxt(filename, self.data.spectrum)

    def save_spectrum_to_file(self, spectrum, saveas):
        filename = os.path.join(self.out_path, saveas)
        logging.info('Spectrum saved to %s ' % filename)
        np.savetxt(filename, spectrum)

    def save_TP_profile(self, save_manual=None, FIT_params=None, FIT_params_std=None, save2pdf=False):
        '''
        function saving TP profiles for individual minimisation/sampling modes.
        If save_manual = True, one can provide FIT_params and FIT_params_std to manually
        save TP_profiles.
        '''
        profilename = '_TP_profile_'
        basename = self.out_path 

        P = self.fitting.forwardmodel.atmosphere.P

        out = np.zeros((len(self.fitting.forwardmodel.atmosphere.P), 3))
        out[:,0] = P # pressure
        
        fit_TPparams_bounds = self.fitting.fit_bounds[self.fitting.fit_X_nparams:]

        if save_manual is not None:
            if len(FIT_params) > 0 and len(FIT_params_std) > 0:
                T_mean, T_sigma, P = iterate_TP_profile(FIT_params, FIT_params_std,fit_TPparams_bounds,
                                                         self.fitting.forwardmodel.atmosphere.TP_profile)

                out[:,1] = T_mean; out[:,2] = T_sigma

                filename = str(basename)+'MANUAL_'+profilename+'.dat'
                logging.info('Saving manual TP_profile to %s' % filename)
                np.savetxt(filename,out)

                if plot:
                    plot_TP_profile(P, T_mean, T_sigma, name='NEST',
                                          save2pdf=save2pdf, out_path=self.out_path)

            else:
                logging.error('Saving manual TP_profile with wrong input parameters')
                exit()

        if self.DOWN:

            fit_params = self.DOWN_params_values
            self.fitting.update_atmospheric_parameters(fit_params)
            T_mean = self.fitting.forwardmodel.T
            T_sigma  = np.zeros_like(T_mean)
            out[:,1] = T_mean
            out[:,2] = T_sigma
            filename = str(basename)+'/DOWN'+profilename+'.dat'
            logging.info('Saving MLE TP_profile to %s' % filename)
            np.savetxt(filename,out)

            if save2pdf:
                plot_TP_profile(P, T_mean, name='DOWN', save2pdf=save2pdf, out_path=self.out_path)
                
            logging.info('Saving MLE contribution function to %s' % (filename))
        
            fit_params = [self.DOWN_out['fit_params'][param]['value'] for param in self.params_names]
            filename2  = str(basename)+'/DOWN_cont_func.dat'
            self.save_contribution_function(fit_params, filename2)

        if self.MCMC:

            logging.info('There are %i different MCMC chains. '
                         'Saving one TP profile for each chain' % len(self.MCMC_out))


            for idx, solution in enumerate(self.MCMC_out):

                #iterate through all upper/lower bounds of parameters to find function minimum and maximum
#                 fit_params = [self.MCMC_out[idx]['fit_params'][param]['value'] for param in self.params_names]

                tp_profile_type = self.atmosphere.TP_type
                if tp_profile_type is 'hybrid' or tp_profile_type is 'rodgers' or tp_profile_type is 'isothermal':
                    tp_iterate = False
                else:
                    tp_iterate = True

                T_mean, T_sigma = iterate_TP_profile(self.MCMC_TP_params_values[idx], self.MCMC_TP_params_std[idx],
                                                     fit_TPparams_bounds, self.fitting.forwardmodel.atmosphere.TP_profile,iterate=tp_iterate)

                out[:,1] = T_mean;
                out[:,2] = T_sigma
                filename = str(basename)+'/MCMC'+profilename+'%i.dat' % (idx)
                logging.info('Saving MCMC TP_profile to %s' % filename)
                np.savetxt(filename,out)

                if save2pdf:
                    plot_TP_profile(P, T_mean, T_sigma, name='MCMC_'+str(idx),
                                          save2pdf=save2pdf, out_path=self.out_path)
                    
                logging.info('Saving MCMC contribution function for solution %i to %s' % (idx, filename))
        
                fit_params = [solution['fit_params'][param]['value'] for param in self.params_names]
                filename2  = str(basename)+'/MCMC_cont_func_%i.dat' % (idx)
                self.save_contribution_function(fit_params, filename2)
                

        if self.NEST:

            logging.info('MultiNest detected %i different modes. '
                         'Saving one TP profile for each solution' % len(self.NEST_out))

            for idx, solution in enumerate(self.NEST_out):

                tp_profile_type = self.atmosphere.TP_type
                if tp_profile_type is 'hybrid' or tp_profile_type is 'rodgers' or tp_profile_type is 'isothermal':
                    tp_iterate = False
                else:
                    tp_iterate = True
                    
                T_mean, T_sigma = iterate_TP_profile(self.NEST_TP_params_values[idx], self.NEST_TP_params_std[idx],
                                                      fit_TPparams_bounds,self.fitting.forwardmodel.atmosphere.TP_profile,iterate=tp_iterate)
            
                out[:,1] = T_mean;
                out[:,2] = T_sigma
                filename = str(basename)+'/NEST'+profilename+'%i.dat' % (idx)

                logging.info('Saving Nested Sampling TP profile for solution %i to %s' % (idx, filename))
                np.savetxt(filename, out)              

                if save2pdf:
                    plot_TP_profile(P, T_mean, T_sigma, name='NEST_'+str(idx),
                                          save2pdf=save2pdf, out_path=self.out_path)

                logging.info('Saving Nested Sampling contribution function for solution %i to %s' % (idx, filename))
        
                fit_params = [solution['fit_params'][param]['value'] for param in self.params_names]
                filename2  = str(basename)+'/NEST_cont_func_%i.dat' % (idx)
                self.save_contribution_function(fit_params, filename2)
                


    def save_contribution_function(self, fit_params,filename):
        
        self.fitting.update_atmospheric_parameters(fit_params)
        tau, tau_total, dtau = self.fitting.forwardmodel.get_contribution_function()
        np.savetxt(filename,tau_total)

    #
    #
    # def analyse_traces(self, tracedata, clr_inv=False, cluster_analysis=False):
    #
    #     ntraces = len(tracedata[:,0])
    #     ncol = len(tracedata[0,:])
    #
    #     if clr_inv == True:
    #         # create new array of traces with the mixing ratios obtained from CLR traces
    #         # log-ratios are always the first n columns
    #         mixing_ratios_clr = tracedata[:, :self.fitting.forwardmodel.atmosphere.nallgases-1] # get traces of log-ratios
    #         mixing_ratios_clr_inv, coupled_mu_trace = self.inverse_clr_transform(mixing_ratios_clr)
    #
    #         self.tracedata_clr_inv = np.c_[mixing_ratios_clr_inv, coupled_mu_trace,
    #                                        tracedata[:, self.fitting.forwardmodel.atmosphere.nallgases-1:]]
    #         tracedata_analyse = self.tracedata_clr_inv # set the trace to analyse and get solutions from
    #
    #         # set parameter names of new traces
    #         self.params_names = []
    #         for idx, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases +
    #                 self.fitting.forwardmodel.atmosphere.inactive_gases):
    #             self.params_names.append(gasname)
    #         ## todo careful about adding coupled_mu, it might interfere with fit-params?
    #         self.params_names.append('coupled_mu')
    #         for i in range(self.fitting.forwardmodel.atmosphere.nallgases-1, len(self.fitting.fit_params_names)):
    #             self.params_names.append(self.fitting.fit_params_names[i])
    #
    #     else:
    #         self.params_names = self.fitting.fit_params_names
    #         tracedata_analyse = tracedata # set the trace to analyse and get solutions from
    #
    #
    #     # @todo: if cluster_analysis is False, but multinest is in multimode, then we should get solutions from multinest output. See end of file with commented code.
    #     if cluster_analysis == True:
    #
    #         # clustering analysis on traces. Split up the traces into individual clusters
    #         logging.info('Cluster analysis with sklearn MeanShift')
    #
    #         from sklearn.cluster import MeanShift
    #         from sklearn.preprocessing import normalize
    #         data_normalized = normalize(tracedata_analyse, axis=0)  # normalize the data
    #         estimator = MeanShift()
    #         estimator.fit(data_normalized)
    #         labels = estimator.labels_
    #
    #         # save txt with labels corresponding to different clusters found in the posterior traces
    #         np.savetxt('Output/labels.txt', labels)
    #         # labels = np.loadtxt('Output/labels.txt')
    #         unique_labels = set(labels)
    #         n_clusters = len(unique_labels) - (1 if -1 in labels else 0)
    #         logging.info('Estimated number of clusters: %d' % n_clusters)
    #
    #     else:
    #         # clustering analysis is not performed. Only one mode assumed
    #         unique_labels = [0]
    #         labels = np.zeros(ntraces)
    #
    #     fit_out = self.get_solutions_from_traces(tracedata_analyse, multimode=cluster_analysis, labels=labels) # todo: needs to be changed, if we get solutions from multinest multimode (see above)
    #     fit_out = self.add_spectra_from_solutions(fit_out)
    #
    #     return fit_out, tracedata_analyse, labels
    #
    # def get_solutions_from_traces(self, tracedata, multimode=False, labels=[]):
    #
    #     ntraces = len(tracedata[:,0])
    #     ncol = len(tracedata[0,:])
    #
    #     if multimode == False:
    #         unique_labels = [0]
    #         labels = np.zeros(ntraces)
    #     else:
    #         unique_labels = set(labels)
    #         if len(labels) != ntraces:
    #             logging.error('Multimode is ON and the number of labels is different from the number of samples')
    #
    #     # Create list of solutions. Each solution (cluster) is stored into a dictionary. If multimode = False, then
    #     # only only one solution is stored.
    #
    #     fit_out = []
    #     for k in unique_labels:
    #         if k >= 0:
    #
    #             if (k > 0 and len(tracedata[labels == k-1, 0])/len(tracedata[labels == k, 0]) > 20):
    #                 # exlcude solutions that have a factor of 20 fewer samples than the previous solution @todo careful!
    #                 break
    #             else:
    #
    #                 # save individual solutions into a dictionary
    #                 dict = {'fit_params': {}}
    #                 sort_order = 0
    #                 for idx, param in enumerate(self.params_names):
    #
    #                     try:
    #                         cluster = tracedata[labels == k, idx] # get cluster of points corresponding to cluster k, and param idx
    #
    #                         # get errors
    #                         hist, bin_edges = np.histogram(cluster, bins=100, density=True)
    #                         bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    #                         idxmax = np.argmax(hist)
    #                         value = bin_centres[idxmax]
    #                         left = hist[:idxmax][::-1]
    #                         right = hist[idxmax:]
    #                         left_centres = bin_centres[:idxmax][::-1]
    #                         right_centres = bin_centres[idxmax:]
    #
    #                         try:
    #                             left_err = value - left_centres[(np.abs(left-np.percentile(left,68))).argmin()]
    #                         except:
    #                             left_err = 0
    #                         try:
    #                             right_err = right_centres[(np.abs(right-np.percentile(right,68))).argmin()] - value
    #                         except:
    #                             right_err = 0
    #                         if left_err == 0 and right_err > 0:
    #                             left_err = right_err
    #                         if right_err == 0 and left_err > 0:
    #                             right_err = left_err
    #
    #                         mean_err = np.mean((left_err, right_err))
    #
    #                         dict['fit_params'][self.params_names[idx]] = {
    #                             'value': value,
    #                             'std': np.std(cluster),
    #                             #'std': mean_err, @ todo try to use mean err!
    #                             'std_plus': right_err,
    #                             'std_minus': left_err,
    #                             'trace': np.asarray(cluster),
    #                             'sort_order': sort_order
    #                         }
    #                         sort_order += 1
    #
    #                     except ValueError:
    #                         # this parameter was not fitted...
    #                         pass
    #                 fit_out.append(dict)
    #
    #     logging.info('Number of solutions identified: %d' % len(fit_out))
    #
    #     return fit_out
    #
