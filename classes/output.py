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
from scipy.optimize import curve_fit
import pickle

#conversion constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg


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

        self.save_fit_output_files()

    #class methods

    def store_downhill_solution(self):

        logging.info('Store the DOWNHILL results')

        self.fitting.update_atmospheric_parameters(self.fitting.DOWN_fit_output)

        dict = {}
        if 'P0' in self.fitting.fit_params_names:
            dict['P0'] = {'value': np.log10(self.fitting.forwardmodel.atmosphere.max_pressure) } # pressure in log space
        if 'Temperature' in self.fitting.fit_params_names:
            dict['Temperature'] = {'value': self.fitting.forwardmodel.atmosphere.planet_temp }
        if 'Radius' in self.fitting.fit_params_names:
            dict['Radius'] = {'value': self.fitting.forwardmodel.atmosphere.planet_radius/RJUP } # pressure in log space
        for idx_gas, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
            dict[gasname] = {'value': self.fitting.forwardmodel.atmosphere.absorbing_gases_X[idx_gas] }
        if not self.params.fit_fix_inactive:
            for idx_gas, gasname in enumerate(self.fitting.forwardmodel.atmosphere.inactive_gases):
                dict[gasname] = {'value': self.fitting.forwardmodel.atmosphere.inactive_gases_X[idx_gas] }

        DOWN_out = [{'fit_params': dict}]

        self.DOWN_out = self.add_spectra_from_solutions(DOWN_out)
        self.DOWN_params_values = self.fitting.DOWN_fit_output
        self.DOWN_TP_params_values = self.DOWN_params_values[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fit_TP_nparams]

        self.params_names = self.fitting.fit_params_names

    def store_mcmc_solutions(self):

        logging.info('Store the MCMC results')

        # read individual traces from PyMCMC output and combine them into a single array, called tracedata
        for thread in glob.glob(os.path.join(self.fitting.dir_mcmc, '*')):
            chainlist = glob.glob(os.path.join(thread, 'Chain_0', 'PFIT*'))
            tracedata = [np.loadtxt(trace) for trace in chainlist]
        tracedata = np.transpose(tracedata)

        if self.params.fit_clr_trans:
            # the mixing ratios are sampled in the log-ratio space. Need to convert them back to simplex
            self.MCMC_out, self.MCMC_tracedata_clr_inv, self.MCMC_labels =  self.analyse_traces(tracedata, multimode=False, clr_inv=True)
            self.MCMC_tracedata = tracedata
        else:
            self.MCMC_out, self.MCMC_tracedata, self.MCMC_labels = self.analyse_traces(tracedata, multimode=False)

        self.MCMC_params_values = [self.MCMC_out[0]['fit_params'][param]['value'] for param in self.params_names]
        self.MCMC_params_std = [self.MCMC_out[0]['fit_params'][param]['std'] for param in self.params_names]
        self.MCMC_X_params_values = self.MCMC_params_values[:self.fitting.fit_X_nparams]
        self.MCMC_X_params_std = self.MCMC_params_std[:self.fitting.fit_X_nparams]
        self.MCMC_TP_params_values = self.MCMC_params_values[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fit_TP_nparams]
        self.MCMC_TP_params_std = self.MCMC_params_std[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fit_TP_nparams]

    def store_nest_solutions(self):

        logging.info('Store the MULTINEST results')

        # get output traces from 1-.txt. Exclude first two columns.
        NEST_data = np.loadtxt(os.path.join(self.fitting.dir_multinest, '1-.txt'))

        self.NEST_tracedata = NEST_data[:,2:] # exclude first two columns
        self.NEST_likelihood = NEST_data[:,:2] # get first two columns

        if self.params.fit_clr_trans:

            # the mixing ratios are sampled in the log-ratio space. Need to convert them back to simplex
            self.NEST_out, self.NEST_tracedata_clr_inv, self.NEST_labels =  self.analyse_traces(self.NEST_tracedata, multimode=True, clr_inv=True)

            # save file with clr traces converted back to simplex, in multinest format
            self.NEST_data_clr_inv = np.c_[NEST_data[:, :2], self.NEST_tracedata_clr_inv]

        else:
            self.NEST_out, self.NEST_tracedata, self.NEST_labels =  self.analyse_traces(self.NEST_tracedata, multimode=True)

        self.NEST_params_values = [self.NEST_out[0]['fit_params'][param]['value'] for param in self.params_names]
        self.NEST_params_std = [self.NEST_out[0]['fit_params'][param]['std'] for param in self.params_names]
        self.NEST_X_params_values = self.NEST_params_values[:self.fitting.fit_X_nparams]
        self.NEST_X_params_std = self.NEST_params_std[:self.fitting.fit_X_nparams]
        self.NEST_TP_params_values = self.NEST_params_values[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fit_TP_nparams]
        self.NEST_TP_params_std = self.NEST_params_std[self.fitting.fit_X_nparams:self.fitting.fit_X_nparams+self.fit_TP_nparams]

    def analyse_traces(self, tracedata, clr_inv=False, multimode=False):

        ntraces = len(tracedata[:,0])
        ncol = len(tracedata[0,:])

        if clr_inv == True:
            # create new array of traces with the mixing ratios obtained from CLR traces
            # log-ratios are always the first n columns
            mixing_ratios_clr = tracedata[:, :self.fitting.forwardmodel.atmosphere.nallgases-1] # get traces of log-ratios
            mixing_ratios_clr_inv, coupled_mu_trace = self.inverse_clr_transform(mixing_ratios_clr)

            self.tracedata_clr_inv = np.c_[mixing_ratios_clr_inv, coupled_mu_trace,
                                           tracedata[:, self.fitting.forwardmodel.atmosphere.nallgases-1:]]
            tracedata_analyse = self.tracedata_clr_inv # set the trace to analyse and get solutions from

            # set parameter names of new traces
            self.params_names = []
            for idx, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases +
                    self.fitting.forwardmodel.atmosphere.inactive_gases):
                self.params_names.append(gasname)
            self.params_names.append('coupled_mu')
            for i in range(self.fitting.forwardmodel.atmosphere.nallgases-1, len(self.fitting.fit_params_names)):
                self.params_names.append(self.fitting.fit_params_names[i])

        else:
            self.params_names = self.fitting.fit_params_names
            tracedata_analyse = tracedata # set the trace to analyse and get solutions from

        if multimode == True:

            # clustering analysis on traces. Split up the traces into individual clusters
            logging.info('Cluster analysis with sklearn MeanShift')

            from sklearn.cluster import MeanShift
            from sklearn.preprocessing import normalize
            data_normalized = normalize(tracedata_analyse, axis=0)  # normalize the data
            estimator = MeanShift()
            estimator.fit(data_normalized)
            labels = estimator.labels_

            # save txt with labels corresponding to different clusters found in the posterior traces
            np.savetxt('Output/labels.txt', labels)
            # labels = np.loadtxt('Output/labels.txt')
            unique_labels = set(labels)
            n_clusters = len(unique_labels) - (1 if -1 in labels else 0)
            logging.info('Estimated number of clusters: %d' % n_clusters)

        else:
            # clustering analysis is not performed. Only one mode assumed
            unique_labels = [0]
            labels = np.zeros(ntraces)

        fit_out = self.get_solutions_from_traces(tracedata_analyse, multimode, labels)
        fit_out = self.add_spectra_from_solutions(fit_out)

        return fit_out, tracedata_analyse, labels

    def get_solutions_from_traces(self, tracedata, multimode=False, labels=[]):

        ntraces = len(tracedata[:,0])
        ncol = len(tracedata[0,:])

        if multimode == False:
            unique_labels = [0]
            labels = np.zeros(ntraces)
        else:
            unique_labels = set(labels)
            if len(labels) != ntraces:
                logging.error('Multimode is ON and the number of labels is different by the number of samples')

        # Create list of solutions. Each solution (cluster) is stored into a dictionary. If multimode = False, then
        # only only one solution is stored.

        # params = ['coupled_mu', 'Temperature', 'Radius', 'P0']
        # params += self.fitting.forwardmodel.atmosphere.absorbing_gases # list of absorbing gases
        # params +=  self.fitting.forwardmodel.atmosphere.inactive_gases # list of inactive gases

        fit_out = []
        for k in unique_labels:
            if k >= 0:

                if (k > 0 and len(tracedata[labels == k-1, 0])/len(tracedata[labels == k, 0]) > 20):
                    # exlcude solutions that have a factor of 20 fewer samples than the previous solution @todo careful!
                    break
                else:

                    # save individual solutions into a dictionary
                    dict = {'fit_params': {}}
                    sort_order = 0
                    for idx, param in enumerate(self.params_names):

                        try:
                            cluster = tracedata[labels == k, idx] # get cluster of points correposning to cluster k, and param idx

                            if param == 'mu' or param == 'coupled_mu':
                                cluster /= AMU # mu in AMU

                            # get errors

                            hist, bin_edges = np.histogram(cluster, bins=100, density=True)
                            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
                            idxmax = np.argmax(hist)
                            value = bin_centres[idxmax]
                            left = hist[:idxmax][::-1]
                            right = hist[idxmax:]
                            left_centres = bin_centres[:idxmax][::-1]
                            right_centres = bin_centres[idxmax:]

                            try:
                                left_err = value - left_centres[(np.abs(left-np.percentile(left,68))).argmin()]
                            except:
                                left_err = 0
                            try:
                                right_err = right_centres[(np.abs(right-np.percentile(right,68))).argmin()] - value
                            except:
                                right_err = 0
                            if left_err == 0 and right_err > 0:
                                left_err = right_err
                            if right_err == 0 and left_err > 0:
                                right_err = left_err

                            mean_err = np.mean((left_err, right_err))

                            dict['fit_params'][self.params_names[idx]] = {
                                'value': value,
                                'std': np.std(cluster),
                                #'std': mean_err, @ todo try to use mean err!
                                'std_plus': right_err,
                                'std_minus': left_err,
                                'trace': np.asarray(cluster),
                                'sort_order': sort_order
                            }
                            sort_order += 1

                        except ValueError:
                            # this parameter was not fitted...
                            pass
                    fit_out.append(dict)

        logging.info('Number of solutions identified: %d' % len(fit_out))

        return fit_out

    def inverse_clr_transform(self, clr_tracedata):

        # Convert the traces of the log-ratios back to the gas mixing ratios in log space. This is done by calculating the inverse
        # CLR in the tracedata array, line by line. The new traces are then returned.

        logging.info('Inverse log-ratio transformation of the traces')

        ntraces = len(clr_tracedata[:,0])
        ncol = len(clr_tracedata[0,:])

        mixing_ratios_tracedata = exp(np.c_[clr_tracedata, -sum(clr_tracedata, axis=1)]) # add log-ratio (= -sum(other log ratios)
        mixing_ratios_tracedata /= (np.asarray([np.sum(mixing_ratios_tracedata, axis=1),]*(ncol+1)).transpose()) # clousure operation

        # calculate couple mean molecular weight trace as a byproduct
        coupled_mu_trace = np.zeros(ntraces)
        i = 0
        for sample in mixing_ratios_tracedata:
            # get mean molecular weight from gas mixing ratios
            absorbing_gases_X = sample[:len(self.fitting.forwardmodel.atmosphere.absorbing_gases)]
            inactive_gases_X = sample[len(self.fitting.forwardmodel.atmosphere.inactive_gases):]
            coupled_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu(absorbing_gases_X, inactive_gases_X)
            coupled_mu_trace[i] = coupled_mu
            i += 1

        return np.power(10, mixing_ratios_tracedata), coupled_mu_trace

    def add_spectra_from_solutions(self, fitting_out):

        for idx, solution in enumerate(fitting_out):

            fit_params = solution['fit_params']

            for idx_gas, gasname in enumerate(self.fitting.forwardmodel.atmosphere.absorbing_gases):
                self.fitting.forwardmodel.atmosphere.absorbing_gases_X[idx_gas] = np.power(10, fit_params[gasname]['value'])
            if not self.params.fit_fix_inactive:
                for idx_gas, gasname in enumerate(self.fitting.forwardmodel.atmosphere.inactive_gases):
                    self.fitting.forwardmodel.atmosphere.inactive_gases_X[idx_gas] = np.power(10, fit_params[gasname]['value'])

            if 'mu' in fit_params:
                self.fitting.forwardmodel.atmosphere.planet_mu = fit_params['mu']['value']*AMU
            else:
                self.fitting.forwardmodel.atmosphere.planet_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu()

            if 'Radius' in fit_params:
                self.fitting.forwardmodel.atmosphere.planet_radius = fit_params['Radius']['value']*RJUP

            if 'P0' in fit_params:
                self.fitting.forwardmodel.atmosphere.max_pressure = np.power(10, fit_params['P0']['value'])

            if 'Temperature' in fit_params:
                self.fitting.forwardmodel.atmosphere.planet_temp = fit_params['Temperature']['value']

            # update atmospheric parameters (PTA, rho, etc)
            self.fitting.forwardmodel.atmosphere.update_atmosphere()

            model = self.fitting.forwardmodel.model()
            model_binned = [model[self.data.spec_bin_grid_idx == i].mean() for i in xrange(1,self.data.n_spec_bin_grid)]

            out = np.zeros((len(self.data.spectrum[:,0]), 2))
            out[:,0] = self.data.spectrum[:,0]
            out[:,1] = model_binned

            fitting_out[idx]['spectrum'] = out

        return fitting_out

    def save_fit_output_files(self):

        # store parameter names of traces
        f = open(os.path.join(self.params.out_path, 'parameters.txt'),'w')
        for param in self.fitting.fit_params_names:
            f.write('%s\n' % param)
        f.close()

        # if clr transform, store parameter names of transformed traces
        if self.params.fit_clr_trans:

            # If clr inverse, write out parameter names with mixing ratio naems
            f = open(os.path.join(self.params.out_path, 'parameters_clr_inv.txt'),'w')
            for param in self.params_names:
                f.write('%s\n' % param)
            f.close()

        if self.DOWN:
            pickle.dump(self.DOWN_out, open('Output/DOWN_out.db', 'wb'))

        if self.MCMC:
            np.savetxt(os.path.join(self.params.out_path, 'MCMC_tracedata.txt'), self.MCMC_tracedata)
            np.savetxt(os.path.join(self.fitting.dir_multinest, 'MCMC_labels.txt'), self.MCMC_labels)
            if self.params.fit_clr_trans:
                np.savetxt(os.path.join(self.params.out_path, 'MCMC_clr_inv_tracedata.txt'), self.MCMC_data_clr_inv)
            pickle.dump(self.MCMC_out, open('Output/MCMC_out.db', 'wb'))

        if self.NEST:
            np.savetxt(os.path.join(self.params.out_path, 'NEST_tracedata.txt'), self.NEST_tracedata)
            np.savetxt(os.path.join(self.params.out_path, 'NEST_likelihood.txt'), self.NEST_likelihood)
            np.savetxt(os.path.join(self.fitting.dir_multinest, 'NEST_labels.txt'), self.NEST_labels)
            if self.params.fit_clr_trans:
                np.savetxt(os.path.join(self.params.out_path, 'NEST_clr_inv_tracedata.txt'), self.NEST_data_clr_inv)
            pickle.dump(self.NEST_out, open('Output/NEST_out.db', 'wb'))

    def plot_all(self, save2pdf=False):

        logging.info('Plotting absolutely everything')

        self.plot_spectrum(save2pdf=save2pdf)
        self.plot_fit(save2pdf=save2pdf)
        self.plot_distributions(save2pdf=save2pdf)

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
        py.plot(self.data.spectrum[:,0], # todo why again?
                self.data.spectrum[:,1],
                color=[0.3,0.3,0.3],
                linewidth=linewidth,label='DATA')

        # plot models
        if self.fitting.DOWN:
            py.plot(self.DOWN_out[0]['spectrum'][:,0], self.DOWN_out[0]['spectrum'][:,1],
                    c='b',label='DOWNHILL',linewidth=linewidth)
        if self.fitting.MCMC:
            for idx, solution in enumerate(self.MCMC_out):
                py.plot(solution['spectrum'][:,0], solution['spectrum'][:,1],
                        label='MCMC %i' % idx, linewidth=linewidth)
        if self.fitting.NEST:
            for idx, solution in enumerate(self.NEST_out):
                py.plot(solution['spectrum'][:,0], solution['spectrum'][:,1],
                        label='NESTED %i' % idx, linewidth=linewidth)

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

    def plot_distributions(self, save2pdf=False):

        logging.info('Plotting sampling distributions. Saving to %s' % self.params.out_path)

        if self.fitting.MCMC:
            plot_posteriors(self.MCMC_out,
                            save2pdf=save2pdf,
                            out_path=self.params.out_path,
                            plot_name = 'mcmc',
                            plot_contour=self.params.out_plot_contour)
        if self.fitting.NEST:
            plot_posteriors(self.NEST_out,
                            save2pdf=save2pdf,
                            out_path=self.params.out_path,
                            plot_name = 'nest',
                            plot_contour=self.params.out_plot_contour)

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

    def save_ascii_spectra(self):

        logging.info('Dumps all model fits and observed spectra to file')

        if self.MCMC:
            for idx, solution in enumerate(self.NEST_out):
                filename = os.path.join(self.params.out_path, 'MCMC_%s_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                logging.info('Saving MCMC spectrum for solution %i to %s' % (idx, filename))
                np.savetxt(filename, solution['spectrum'])

        if self.NEST:
            logging.info('MultiNest detected %i different modes. '
                         'Saving one model spectrum for each solution' % len(self.NEST_out))
            for idx, solution in enumerate(self.NEST_out):
                filename = os.path.join(self.params.out_path, 'NEST_%s_spectrum_%i.dat' % (self.__MODEL_ID__, idx))
                logging.info('Saving Nested Sampling spectrum for solution %i to %s' % (idx, filename))
                np.savetxt(filename, solution['spectrum'])
        if self.DOWN:
                filename = os.path.join(self.params.out_path, 'DOWN_%s_spectrum.dat' % (self.__MODEL_ID__))
                logging.info('Saving DOWHNILL spectrum to %s' % filename)
                np.savetxt(filename, self.DOWN_out[0]['spectrum'])

        filename = os.path.join(self.params.out_path, 'observed_%s_spectrum.dat' % (self.__MODEL_ID__))
        np.savetxt(filename, self.data.spectrum)


    # the mixing ratio are not sampled in the log-ratio space, hence we can use the different
        # modes detected by multinest and get individual solutions, errors and chains
        #     self.params_names = self.fitting.fit_params_names
        #
        #     NEST_analyzer = pymultinest.Analyzer(n_params=len(self.fitting.fit_params),
        #                                          outputfiles_basename=os.path.join(self.fitting.dir_multinest, '1-'))
        #     NEST_stats = NEST_analyzer.get_stats()
        #     self.NEST_out = []
        #
        #     modes = []
        #     chains = []
        #
        #     if self.params.nest_multimodes:
        #         # get separate chains for each mode
        #         with open(os.path.join(self.fitting.dir_multinest, '1-post_separate.dat')) as f:
        #             lines = f.readlines()
        #             for idx, line in enumerate(lines):
        #                 if idx > 2: # skip the first two lines
        #                     if lines[idx-1] == '\n' and lines[idx-2] == '\n':
        #                         modes.append(chains)
        #                         chains = []
        #                 chain = [float(x) for x in line.split()[2:]]
        #
        #                 if len(chain) > 0:
        #                     chains.append(chain)
        #             modes.append(chains)
        #         modes_array = []
        #         for mode in modes:
        #             mode_array = np.zeros((len(mode), len(mode[0])))
        #             for idx, line in enumerate(mode):
        #                 mode_array[idx, :] = line
        #             modes_array.append(mode_array)
        #     else:
        #         # not running in multimode. Get chain directly from file 1-.txt
        #         data = np.loadtxt(os.path.join(self.fitting.dir_multinest, '1-.txt'))
        #         modes_array = [data[:,2:]]
        #
        #     import pickle
        #     pickle.dump(modes_array, open('Output/modes.db', 'wb'))
        #
        #     for n in range(len(modes)):
        #
        #         dict = {'fit_params': {}}
        #         sort_order = 0
        #         for i in range(len(self.fitting.fit_params)):
        #             value = NEST_stats['modes'][n]['maximum a posterior'][i]
        #             std = NEST_stats['modes'][n]['sigma'][i]
        #             chain = modes_array[n][:,i]
        #
        #             if self.params_names[i] == 'mu':
        #                 value /= AMU
        #                 std /= AMU
        #                 chain /= AMU
        #             elif self.params_names[i] == 'P0':
        #                 value = np.power(10, mu)
        #                 std *= value
        #                 chain = np.power(10, chain)
        #
        #             dict['fit_params'][self.params_names[i]] = {
        #                 'value': value,
        #                 'std': std,
        #                 'chain': chain,
        #                 'sort_order': sort_order
        #             }
        #             sort_order += 1
        #
        #         self.NEST_out.append(dict)
        #
        # self.NEST_out = self.add_spectra_from_solutions(self.NEST_out)
        # pickle.dump(self.MCMC_out, open('Output/NEST_out.db', 'wb'))
