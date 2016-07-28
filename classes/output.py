'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Output class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

# loading libraries
import copy
import ctypes as C
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle

from library_constants import *
from library_general import *
import logging

import matplotlib.pylab as plt

try:
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

class output(object):

    def __init__(self, fitting=None, forwardmodel=None, data=None, atmosphere=None, params=None, out_path=None):

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

        if fitting is not None:

            types = ['downhill', 'mcmc', 'nest']
            func_types = [self.store_downhill_solution,
                          self.store_mcmc_solutions,
                          self.store_nest_solutions]
            run = [self.fitting.DOWN, self.fitting.MCMC, self.fitting.NEST]
            out_filenames = [self.params.downhill_out_filename,
                             self.params.mcmc_out_filename,
                             self.params.nest_out_filename]
            self.dbfilename = {}
            for idx, val in enumerate(types):
                if run[idx]:

                    logging.info('Storing %s solution' % val)

                    outdb = self.initialize_output(val)
                    outdb = func_types[idx](outdb)
                    outdb = self.add_data_from_solutions(outdb)

                    if out_filenames[idx] == 'default':
                        filename = os.path.join(self.params.out_path, '%s_out.db' % types[idx])
                    else:
                        filename = out_filenames[idx]
                    self.dbfilename[val] = filename
                    pickle.dump(outdb, open(filename, 'wb'), protocol=2)

        logging.info('Output object correctly initialised')

    def initialize_output(self, fitting_type):

        outdb = {'type': fitting_type,
                 'obs_spectrum': self.data.obs_spectrum,
                 'params': self.params.params_to_dict(),
                 'fit_params_names': self.fitting.fit_params_names,
                 'fit_params_bounds': self.fitting.fit_bounds,
                 'fit_params_texlabels': self.fitting.fit_params_texlabels,
                 'solutions': []}

        if os.path.isfile(self.params.in_spectrum_db):


            try:
                SPECTRUM_db = pickle.load(open(self.params.in_spectrum_db, 'rb'), encoding='latin1')
            except:
                SPECTRUM_db = pickle.load(open(self.params.in_spectrum_db))
            outdb['SPECTRUM_db'] = SPECTRUM_db

        return outdb

    def store_downhill_solution(self, DOWN_out):

        logging.info('Store the downhill results')

        dict = {}
        for idx, param_name in enumerate(self.fitting.fit_params_names):
            dict[param_name] = {'value': self.fitting.DOWN_fit_output[idx] }

        DOWN_out['solutions'].append({'type': 'down', 'fit_params': dict})

        return DOWN_out

    def store_mcmc_solutions(self, MCMC_out):

        logging.info('Store the MCMC results')

        # read individual traces from PyMCMC output and combine them into a single array
        for thread in glob.glob(os.path.join(self.fitting.dir_mcmc, '*')):
            chainlist = glob.glob(os.path.join(thread, 'Chain_0', 'PFIT*'))
            tracedata = [np.loadtxt(trace) for trace in chainlist]
        tracedata = np.transpose(tracedata)

        # todo MCMC output not supported at the moment...

    def store_nest_solutions(self, NEST_out):

        logging.info('Store the multinest results')

        data = np.loadtxt(os.path.join(self.fitting.dir_multinest, '1-.txt'))

        NEST_analyzer = pymultinest.Analyzer(n_params=len(self.fitting.fit_params),
                                             outputfiles_basename=os.path.join(self.fitting.dir_multinest, '1-'))
        NEST_stats = NEST_analyzer.get_stats()
        NEST_out['NEST_stats'] = NEST_stats
        NEST_out['global_logE'] = (NEST_out['NEST_stats']['global evidence'], NEST_out['NEST_stats']['global evidence error'])

        modes = []
        modes_weights = []
        chains = []
        chains_weights = []

        if self.params.nest_multimodes:

            # separate modes. get individual samples for each mode

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

        modes_weights = np.asarray(modes_weights)

        for nmode in range(len(modes)):

            dict = {'type': 'nest',
                    'local_logE': (NEST_out['NEST_stats']['modes'][0]['local log-evidence'],  NEST_out['NEST_stats']['modes'][0]['local log-evidence error']),
                    'weights': modes_weights[nmode],
                    'tracedata': modes_array[nmode],
                    'fit_params': {}}

            for idx, param_name in enumerate(self.fitting.fit_params_names):

                trace = modes_array[nmode][:,idx]
                q_16, q_50, q_84 = quantile_corner(trace, [0.16, 0.5, 0.84],
                            weights=np.asarray(modes_weights[nmode]))
                dict['fit_params'][param_name] = {
                    'value' : q_50,
                    'sigma_m' : q_50-q_16,
                    'sigma_p' : q_84-q_50,
                    'nest_map': NEST_stats['modes'][nmode]['maximum a posterior'][idx],
                    'nest_mean': NEST_stats['modes'][nmode]['mean'][idx],
                    'nest_sigma': NEST_stats['modes'][nmode]['sigma'][idx],
                    'trace': trace,
                }

            if self.params.fit_clr_trans == True:

                # todo check!

                nallgases = self.fitting.forwardmodel.atmosphere.nallgases

                # if centered-log-ratio transformation is True, transform the traces of the log ratios back to abundances
                mixing_ratios_clr = modes_array[nmode][:, :nallgases-1] # get traces of log-ratios
                mixing_ratios_clr = np.c_[mixing_ratios_clr, -sum(mixing_ratios_clr, axis=1)] #add last clr (= -sum(clr) )
                mixing_ratios_clr_inv, coupled_mu_trace = inverse_clr_transform(mixing_ratios_clr)

                self.NEST_data_clr_inv = np.c_[mixing_ratios_clr_inv, coupled_mu_trace]

                # get means
                clr = np.asarray([NEST_stats['modes'][nmode]['maximum a posterior'][i] for i in range(nallgases-1)])
                clr = np.append(clr, -np.sum(clr))

                # add last CLR to dictionary. This is not fitted, as it is equal to minus the sum of all the other CLRs
                # Note that that the sum of all CLR is equal to zero
                dict['fit_params']['%s_CLR' % self.forwardmodel.atmosphere.inactive_gases[-1]] = {
                    'map': clr[nallgases-1],
                    'sigma': weighted_avg_and_std(mixing_ratios_clr[:, nallgases-1], modes_weights[nmode])[1], # weighted std
                    'trace': mixing_ratios_clr[:, nallgases-1],
                }

                mixing_means = np.exp(clr) # add log-ratio (= -sum(other log ratios)
                mixing_means /= sum(mixing_means) # closure operation

                # get coupled std of means
                absorbing_gases_X = mixing_means[:len(self.params.atm_active_gases)]
                inactive_gases_X = mixing_means[len(self.params.atm_inactive_gases)+1:]

                coupled_mu = self.fitting.forwardmodel.atmosphere.get_coupled_planet_mu(absorbing_gases_X, inactive_gases_X)

                # convert mixing means to log space # todo consider new parameter X_log
                mixing_means = np.log10(mixing_means)

                # set new list of params names, including the clr inv mixing ratios
                self.clrinv_params_names = []
                for idx, gasname in enumerate(self.params.atm_active_gases +
                        self.params.atm_inactive_gases):
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
                        'sigma': std,
                        'trace': trace,
                    }

            NEST_out['solutions'].append(dict)

        return NEST_out

    def add_data_from_solutions(self, fitting_out):


        for idx, solution in enumerate(fitting_out['solutions']):

            logging.info('Solution %i' % idx)

            fit_params = [solution['fit_params'][param]['value'] for param in self.fitting.fit_params_names]
            solution = fitting_out['solutions'][idx]

            solution = self.get_spectra(solution)  # compute spectra, contribution from opacities, and contribution function
            solution = self.get_profiles(solution) # compute mixing ratio and tp profiles

            fitting_out['solutions'][idx] = solution

        return fitting_out

    def get_spectra(self, solution):

        fit_params = [solution['fit_params'][param]['value'] for param in self.fitting.fit_params_names]

        # load model using full wavenumber range
        if self.params.gen_type is 'transmission' or self.params.fit_transmission:
            self.fitting.forwardmodel.atmosphere.load_opacity_arrays(wngrid='full')
            wavegrid  = self.data.int_wlgrid_full
            nwavegrid = self.data.int_nwlgrid_full
            bingrid   = self.data.intsp_bingrididx_full
            nbingrid = self.data.intsp_nbingrid_full
        else:
            self.fitting.forwardmodel.atmosphere.load_opacity_arrays(wngrid='obs_spectrum')
            wavegrid  = self.data.int_wlgrid_obs
            nwavegrid = self.data.int_nwlgrid_obs
            bingrid   = self.data.intsp_bingrididx
            nbingrid = self.data.intsp_nbingrid

        # update atmospheric parameters to current solution
        self.fitting.update_atmospheric_parameters(fit_params)

        model = self.fitting.forwardmodel.model()
        model[np.isnan(model)] = 0.0
        model[np.isinf(model)] = 0.0


        # store observed spectrum
        solution['obs_spectrum'] = np.zeros((self.data.obs_nwlgrid, 6))
        solution['obs_spectrum'][:,0] = self.data.obs_spectrum[:,0]
        solution['obs_spectrum'][:,1] = self.data.obs_spectrum[:,1]
        solution['obs_spectrum'][:,2] = self.data.obs_spectrum[:,2]

        #append fitted spectrum to observed spectrum array
        solution['obs_spectrum'][:,3] = np.asarray([model[bingrid == i].mean()
                                                   for i in range(1, nbingrid+1)])

        solution['fit_spectrum_xsecres'] = np.zeros((nwavegrid, 3))
        solution['fit_spectrum_xsecres'][:,0] = wavegrid
        solution['fit_spectrum_xsecres'][:,1] = model


        # calculate 1 sigma spectrum
        if self.params.out_sigma_spectrum and solution['type'] == 'nest':
            sigmasp = self.get_one_sigma_spectrum(solution)
            solution['fit_spectrum_xsecres'][:,2] = sigmasp
            solution['obs_spectrum'][:,4] = np.asarray([sigmasp[bingrid == i].mean()
                                                        for i in range(1, nbingrid+1)])

        # calculate contribution function
        contrib_func = self.fitting.forwardmodel.model(return_tau=True)
        solution['contrib_func'] = contrib_func

        # calculate spectral contribution from individual opacity sources
        solution['opacity_contrib'] = {}
        active_mixratio_profile = self.fitting.forwardmodel.atmosphere.active_mixratio_profile
        atm_rayleigh = self.fitting.forwardmodel.params.atm_rayleigh
        atm_cia = self.fitting.forwardmodel.params.atm_cia
        atm_clouds = self.fitting.forwardmodel.params.atm_clouds

        # opacity from molecules
        for idx, val in enumerate(self.atmosphere.active_gases):
            mask = np.ones(len(self.atmosphere.active_gases), dtype=bool)
            mask[idx] = 0

            active_mixratio_profile_mask = np.copy(active_mixratio_profile)
            active_mixratio_profile_mask[mask, :] = 0
            self.fitting.forwardmodel.atmosphere.active_mixratio_profile = active_mixratio_profile_mask
            self.fitting.forwardmodel.params.atm_rayleigh = False
            self.fitting.forwardmodel.params.atm_cia = False
            #self.fitting.forwardmodel.params.atm_clouds = False
            solution['opacity_contrib'][val] = np.zeros((nwavegrid, 2))
            solution['opacity_contrib'][val][:,0] = wavegrid
            solution['opacity_contrib'][val][:,1] = self.fitting.forwardmodel.model(mixratio_mask=mask)
            

        self.fitting.forwardmodel.atmosphere.active_mixratio_profile = np.copy(active_mixratio_profile)

        if self.params.gen_type == 'transmission':

            self.fitting.forwardmodel.atmosphere.active_mixratio_profile[:, :] = 0

            # opacity from rayleigh
            if atm_rayleigh:
                self.fitting.forwardmodel.params.atm_rayleigh = True
                self.fitting.forwardmodel.params.atm_cia = False
                #self.fitting.forwardmodel.params.atm_clouds = False
                solution['opacity_contrib']['rayleigh'] = np.zeros((self.data.int_nwlgrid_full, 2))
                solution['opacity_contrib']['rayleigh'][:,0] = self.data.int_wlgrid_full
                solution['opacity_contrib']['rayleigh'][:,1] = self.fitting.forwardmodel.model()

            # opacity from cia
            if atm_cia:
                self.fitting.forwardmodel.params.atm_rayleigh = False
                self.fitting.forwardmodel.params.atm_cia = True
                #self.fitting.forwardmodel.params.atm_clouds = False
                solution['opacity_contrib']['cia'] = np.zeros((self.data.int_nwlgrid_full, 2))
                solution['opacity_contrib']['cia'][:,0] = self.data.int_wlgrid_full
                solution['opacity_contrib']['cia'][:,1] = self.fitting.forwardmodel.model()

            # opacity from clouds
            if atm_clouds:
                self.fitting.forwardmodel.params.atm_rayleigh = False
                self.fitting.forwardmodel.params.atm_cia = False
                self.fitting.forwardmodel.params.atm_clouds = True
                solution['opacity_contrib']['clouds'] = np.zeros((self.data.int_nwlgrid_full, 2))
                solution['opacity_contrib']['clouds'][:,0] = self.data.int_wlgrid_full
                solution['opacity_contrib']['clouds'][:,1] = self.fitting.forwardmodel.model()

            self.fitting.forwardmodel.atmosphere.active_mixratio_profile = np.copy(active_mixratio_profile)

        self.fitting.forwardmodel.params.atm_rayleigh = atm_rayleigh
        self.fitting.forwardmodel.params.atm_cia = atm_cia
        self.fitting.forwardmodel.params.atm_clouds = atm_clouds

        return solution

    def get_profiles(self, solution):

        logging.info('Store TP profile and mixing ratios')

        # best solution
        fit_params = [solution['fit_params'][param]['value'] for param in self.fitting.fit_params_names]

        if solution['type'] == 'nest':
            # compute standard deviation for mixing ratio and tp profiles
            std_tpprofiles, std_molprofiles_active, std_molprofiles_inactive = self.get_one_sigma_profiles(solution)

        # update atmospheric parameters to best solution
        self.fitting.update_atmospheric_parameters(fit_params)

        # tp profile
        solution['tp_profile'] = np.zeros((self.atmosphere.nlayers, 3))
        solution['tp_profile'][:,0] = self.fitting.forwardmodel.atmosphere.pressure_profile
        solution['tp_profile'][:,1] = self.fitting.forwardmodel.atmosphere.temperature_profile
        if solution['type'] == 'nest':
            solution['tp_profile'][:,2] = std_tpprofiles

        # mixing ratios
        solution['active_mixratio_profile'] = np.zeros((len(self.atmosphere.active_gases), self.atmosphere.nlayers, 3))
        solution['inactive_mixratio_profile'] = np.zeros((len(self.atmosphere.inactive_gases), self.atmosphere.nlayers, 3))
        for i in range(len(self.atmosphere.active_gases)):
            solution['active_mixratio_profile'][i,:,0] = self.fitting.forwardmodel.atmosphere.pressure_profile
            solution['active_mixratio_profile'][i,:,1] =  self.fitting.forwardmodel.atmosphere.active_mixratio_profile[i,:]
            if solution['type'] == 'nest':
                solution['active_mixratio_profile'][i,:,2] =  std_molprofiles_active[i,:]
        for i in range(len(self.atmosphere.inactive_gases)):
            solution['inactive_mixratio_profile'][i,:,0] = self.fitting.forwardmodel.atmosphere.pressure_profile
            solution['inactive_mixratio_profile'][i,:,1] =  self.fitting.forwardmodel.atmosphere.inactive_mixratio_profile[i,:]
            if solution['type'] == 'nest':
                solution['inactive_mixratio_profile'][i,:,2] =  std_molprofiles_inactive[i,:]

        return solution

    def get_one_sigma_spectrum(self, solution):

        logging.info('Compute 1 sigma model spectrum spread')

        # best solution
        fit_params = [solution['fit_params'][param]['value'] for param in self.fitting.fit_params_names]

        weights = []
        nspectra = int(self.params.out_sigma_spectrum_frac * np.shape(solution['tracedata'])[0])

        models = np.zeros((nspectra, self.data.int_nwlgrid_full)) # number of possible combinations
        weights = np.zeros((nspectra))
        for i in range(nspectra):
            rand_idx = random.randint(0, np.shape(solution['tracedata'])[0])
            fit_params_iter = solution['tracedata'][rand_idx]
            weights[i] = solution['weights'][i]
            self.fitting.update_atmospheric_parameters(fit_params_iter)
            model = self.fitting.forwardmodel.model()
            models[i, :] = self.fitting.forwardmodel.model()

        models = models[1:,:] # exclude 1st spectrum, for some reasons is made of nan
        weights = np.asarray(weights[1:])

        std_spectrum = np.zeros((self.data.int_nwngrid_full))
        for i in range(self.data.int_nwngrid_full):
            std_spectrum[i] =  weighted_avg_and_std(models[:,i], weights=weights, axis=0)[1]

        # update atmospheric parameters to best solution
        self.fitting.update_atmospheric_parameters(fit_params)

        return std_spectrum


    def get_one_sigma_profiles(self, solution):

        logging.info('Compute 1 sigma mixing ratios and tp profiles')

        sol_tracedata = solution['tracedata']
        sol_weights = solution['weights']
        nprofiles = len(sol_tracedata)
        tpprofiles = np.zeros((nprofiles, self.atmosphere.nlayers))
        molprofiles_active = np.zeros((nprofiles, self.atmosphere.nactivegases, self.atmosphere.nlayers))
        molprofiles_inactive = np.zeros((nprofiles, self.atmosphere.ninactivegases, self.atmosphere.nlayers))
        for i in range(nprofiles):
            fit_params_iter = sol_tracedata[i]
            self.fitting.update_atmospheric_parameters(fit_params_iter)
            tpprofiles[i, :] = self.fitting.forwardmodel.atmosphere.temperature_profile
            for j in range(self.atmosphere.nactivegases):
                molprofiles_active[i,j,:] = self.atmosphere.active_mixratio_profile[j,:]
            for j in range(self.atmosphere.ninactivegases):
                molprofiles_inactive[i,j,:] = self.atmosphere.inactive_mixratio_profile[j,:]
            std_tpprofiles = np.zeros((self.atmosphere.nlayers))
        std_molprofiles_active = np.zeros((self.atmosphere.nactivegases, self.atmosphere.nlayers))
        std_molprofiles_inactive = np.zeros((self.atmosphere.ninactivegases, self.atmosphere.nlayers))
        for i in range(self.atmosphere.nlayers):
            std_tpprofiles[i] = weighted_avg_and_std(tpprofiles[:,i], weights=sol_weights, axis=0)[1]
            for j in range(self.atmosphere.nactivegases):
                std_molprofiles_active[j,i] = weighted_avg_and_std(molprofiles_active[:,j,i], weights=sol_weights, axis=0)[1]
            for j in range(self.atmosphere.ninactivegases):
                std_molprofiles_inactive[j,i] = weighted_avg_and_std(molprofiles_inactive[:,j,i], weights=sol_weights, axis=0)[1]

        return std_tpprofiles, std_molprofiles_active, std_molprofiles_inactive
