'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Atmosphere class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries
import numpy as np
import scipy.special as spe
import logging
import sys
import ctypes as C

from multiprocessing import Process, Queue
import time

import matplotlib.pylab as plt

from library_constants import *
from library_general import *

try:
    import library_cythonised_functions as cy_fun
    cythonised = True
except ImportError:
    cythonised = False

cythonised = False

class atmosphere(object):

    def __init__(self, data, tp_profile_type=None, covariance=None, nthreads=1):

        logging.info('Initialising atmosphere object')

        self.params = data.params

        self.nthreads = nthreads

        self.data = data

        # set planet radius, mass, gravity
        self.planet_radius = self.params.planet_radius
        self.planet_mass = self.params.planet_mass

        self.nlayers = self.params.atm_nlayers
        self.nlevels = self.nlayers + 1

        # set pressure profile of layer boundaries
        self.set_pressure_profile()

        logging.info('Atmospheric pressure boundaries from internal model: %f-%f' %
                     (np.min(self.pressure_profile), np.max(self.pressure_profile)))

        if self.params.atm_active_gases[0] == 'FILE':
            self.set_gases_from_file()
        else:
            self.set_gases_from_param()

        self.nallgases = len(self.active_gases) + len(self.inactive_gases)
        self.nactivegases = len(self.active_gases)
        self.ninactivegases = len(self.inactive_gases)

        # set planet mean molecular weight
        self.set_mu_profile()

        # set the temperature profile
        if self.params.atm_tp_type.upper() == 'FILE':
            self.temperature_profile  = self.get_temperature_profile_from_file()
        else:
            self.temperature_profile  = self.get_temperature_profile_from_param()

        # compute altitude profile, scale height and planet gravity arrays
        self.set_altitude_gravity_scaleheight_profile()

        # set density profile
        self.set_density_profile()

        logging.info('Star radius: %.3f RSUN = %.9e m' % (self.params.star_radius/RSOL, self.params.star_radius))
        logging.info('Planet radius: %.3f RJUP = %.9e m' % (self.planet_radius/RJUP, self.planet_radius))
        logging.info('Planet mass: %.3f MJUP = %.9e kg' % (self.planet_mass/MJUP, self.planet_mass))
        logging.info('Planet gravity (log g) (1st layer): %.3f cgs = %.3f m/s2' % (np.log10(self.gravity_profile[0]*100.),
                                                                                   self.gravity_profile[0]))
        logging.info('Mean molecular weight (1st layer): %.5f AMU' % (self.mu_profile[0]/AMU))
        logging.info('Scale height (1st layer): %.1f km' % (self.scaleheight_profile[0]/1000.))
        logging.info('Temperature (1st layer): %.1f K' % (self.temperature_profile[0]))
        logging.info('Atmospheric max pressure: %.3f bar' % (self.params.atm_max_pres/1e5))

        # selecting TP profile to use
        if tp_profile_type is None:
            self.TP_type = self.params.atm_tp_type
        else:
            self.TP_type = tp_profile_type

        # set non-linear sampling grid for hybrid model
        if tp_profile_type == 'hybrid':
            self.hybrid_covmat = covariance
            self.get_TP_sample_grid(covariance, delta=0.05)

        # load opacity arrays for the appropriate wavenumber grid (gas, rayleigh, cia)
        self.opacity_wngrid = ''
        if self.params.mode == 'retrieval':
            self.load_opacity_arrays(wngrid='obs_spectrum', nthreads=nthreads)
        else:
            if self.params.gen_manual_waverange:
                self.load_opacity_arrays(wngrid='manual', nthreads=nthreads)
            else:
                self.load_opacity_arrays(wngrid='native', nthreads=nthreads)

        # initialise ACE specific parameters
        if self.params.gen_ace:

            # loading Fortran code for chemically consistent model
            self.ace_lib = C.CDLL('./library/ACE/ACE.so', mode=C.RTLD_GLOBAL)

            # set solar elemental abundances. H and He are always set to solar and never change.
            self.ace_H_solar = 12.
            self.ace_He_solar = 10.93
            self.ace_C_solar = 8.43
            self.ace_O_solar = 8.69
            self.ace_N_solar = 7.83
            self.ace_metallicity = self.params.atm_ace_metallicity
            self.ace_co = self.params.atm_ace_co
            self.set_ace_params()

        logging.info('Atmosphere object initialised')

    def set_gases_from_file(self):

        # loading mixing ratio profiles from external file (preloaded in data object)
        self.active_gases = self.data.active_gases_file
        self.inactive_gases = ['H2', 'HE', 'N2']

        # interpolate mixing ratio profiles to pressure grid
        self.active_mixratio_profile = np.zeros((len(self.active_gases), self.nlayers))
        for idx in range(len(self.active_gases)):
            self.active_mixratio_profile[idx,:] = np.interp(self.pressure_profile,
                                                            self.data.mixratio_pressure_file,
                                                            self.data.active_mixratio_profile_file[:,idx])
        self.inactive_mixratio_profile = np.zeros((len(self.inactive_gases), self.nlayers))
        for idx in range(len(self.inactive_gases)):
            self.inactive_mixratio_profile[idx,:] = np.interp(self.pressure_profile,
                                                              self.data.mixratio_pressure_file,
                                                              self.data.inactive_mixratio_profile_file[:,idx])

    def set_gases_from_param(self):

        # active and inactive gases set from param file

        # active and inactive gas names
        self.active_gases = self.params.atm_active_gases
        self.inactive_gases = ['H2', 'HE', 'N2'] # using capital letters! We only consider H2, He and N2.

        # set mixing ratios profiles of active gases
        self.active_mixratio_profile = np.zeros((len(self.active_gases), self.nlayers))
        for i in range(len(self.active_gases)):
            self.active_mixratio_profile[i, :] = self.params.atm_active_gases_mixratios[i]

        # Derive mixing ratio of inactive gases.

        # There are three inactive gases, set in this order in all arrays: H2, HE, N2
        self.inactive_mixratio_profile = np.zeros((len(self.inactive_gases), self.nlayers))

        # set N2 mixing ratio
        self.inactive_mixratio_profile[2, :] = self.params.atm_N2_mixratio

        # set H2 and He mixing ratios

        # first get the sum of the mixing ratio of all active gases
        if len(self.active_gases) > 1:
            active_mixratio_sum = np.sum(self.active_mixratio_profile, axis = 0)
        else:
            active_mixratio_sum = np.copy(self.active_mixratio_profile)


        # add the N2 mixing ratio profile to this sum
        active_mixratio_sum += self.inactive_mixratio_profile[2, :]

        # error if this sum is larger than 1 (i.e. 100%)
        if np.any(active_mixratio_sum > 1.):
            logging.error('It seems that the sum of the mixing ratios of the active gases is larger than unity in at least '
                          'one atmospheric layer. Solve this problem before continuing!')
            exit()

        # The remainder of the atmosphere is made of a mixture of He/H2 , with ratio given by the He_H2_ratio parameter
        mixratio_remainder = 1. - active_mixratio_sum
        self.inactive_mixratio_profile[0, :] = mixratio_remainder/(1. + self.params.atm_He_H2_ratio) # H2
        self.inactive_mixratio_profile[1, :] =  self.params.atm_He_H2_ratio * self.inactive_mixratio_profile[0, :]


    def load_opacity_arrays(self, wngrid='obs_spectrum', nthreads=1):

        if self.opacity_wngrid != wngrid:

            logging.info('Loading opacity arrays for grid `%s`' % wngrid)

            self.opacity_wngrid = wngrid

            # select the wavenumber grid. The grids are defined in the data class, function load_wavenumber_grid
            if wngrid == 'obs_spectrum':
                self.int_nwngrid = self.data.int_nwngrid_obs
                self.int_wngrid = self.data.int_wngrid_obs
                self.int_wngrid_idxmin = self.data.int_wngrid_obs_idxmin
                self.int_wngrid_idxmax = self.data.int_wngrid_obs_idxmax

            elif wngrid == 'native':
                self.int_nwngrid = self.data.int_nwngrid_native
                self.int_wngrid = self.data.int_wngrid_native
                self.int_wngrid_idxmin = 0
                self.int_wngrid_idxmax = self.data.int_nwngrid_native

            elif wngrid == 'manual':
                self.int_nwngrid = self.data.int_nwngrid_manual
                self.int_wngrid = self.data.int_wngrid_manual
                self.int_wngrid_idxmin = self.data.int_wngrid_manual_idxmin
                self.int_wngrid_idxmax = self.data.int_wngrid_manual_idxmax

            elif wngrid == 'extended':
                self.int_nwngrid = self.data.int_nwngrid_extended
                self.int_wngrid = self.data.int_wngrid_extended
                self.int_wngrid_idxmin = self.data.int_wngrid_extended_idxmin
                self.int_wngrid_idxmax = self.data.int_wngrid_extended_idxmax
            else:

                logging.error('Cannot load the opacity arrays for grid `%s`' % wngrid)
                exit()

            self.int_wlgrid = 10000./self.int_wngrid
            self.int_nwlgrid = self.int_nwngrid


            if self.params.mode == 'retrieval':
                # get bins and indexes for spectrum binning (used only if opacity method is xsec, not used for ktables)
                self.int_bingrid, self.int_bingrididx = get_specbingrid(self.data.obs_wlgrid,
                                                                         self.int_wlgrid,
                                                                         self.data.obs_binwidths)
                self.int_nbingrid = len(self.data.obs_binwidths)

            # load SED array for emission
            if self.params.gen_type.upper() == 'EMISSION':
                self.star_sed =  self.data.star_sed_native[self.int_wngrid_idxmin:self.int_wngrid_idxmax]

            # load arrays (interpolate to pressure profile and restrict wavenumber range to selected wngrid)

            if self.params.in_opacity_method in ['xsec_lowres', 'xsec_highres', 'xsec_sampled', 'xsec']:
                # get sigma array (and interpolate sigma array to pressure profile)
                self.sigma_array = self.get_sigma_array(nthreads=nthreads)
                self.sigma_array_flat = self.sigma_array.flatten()
            elif self.params.in_opacity_method in ['ktab', 'ktable', 'ktables']:
                # get sigma array (and interpolate sigma array to pressure profile)
                self.ktables_array = self.get_ktables_array(nthreads=nthreads)
                self.ktables_array_flat = self.ktables_array.flatten()

            # get sigma rayleigh array (for active and inactive absorbers)
            self.sigma_rayleigh_array = self.get_sigma_rayleigh_array()
            self.sigma_rayleigh_array_flat = self.sigma_rayleigh_array.flatten()

            # get collision induced absorption cross sections
            self.sigma_cia_array = self.get_sigma_cia_array()
            self.sigma_cia_array_flat = self.sigma_cia_array.flatten()

            # get the gas indexes of the molecules inside the pairs
            self.cia_idx = self.get_cia_idx()

            # get clouds specific parameters
            self.clouds_pressure = self.params.atm_clouds_pressure
            
            #get mie scattrering specific parameters
            self.mie_r = self.params.atm_mie_r
            self.mie_q = self.params.atm_mie_q
            self.mie_f = self.params.atm_mie_f
            self.mie_topP    = self.params.atm_mie_topP
            self.mie_bottomP = self.params.atm_mie_bottomP
            
            self.sigma_mie_array = self.get_sigma_mie_array()
#             print self.sigma_mie_array
#             plt.figure()
#             plt.plot(self.sigma_mie_array)
#             plt.show()
#             
#             exit()

        else:
            logging.info('Opacity for grid `%s` already loaded' % wngrid)

    def set_mu_profile(self):

        # get mu for each layer
        mu = np.zeros(self.nlayers)
        for i in range(self.nlayers):
            for idx, gasname in enumerate(self.active_gases):
                mu[i] += self.active_mixratio_profile[idx, i] * get_molecular_weight(gasname)
            for idx, gasname in enumerate(self.inactive_gases):
                mu[i] += self.inactive_mixratio_profile[idx, i] * get_molecular_weight(gasname)
            #logging.debug('Mean molecular weight for layer %i is %.4f' % (i, mu[i]/AMU))

        self.mu_profile = mu


    # get the pressure profile
    def set_pressure_profile(self):

        # set pressure profile of layer boundaries
        press_exp = np.linspace(np.log(self.params.atm_min_pres), np.log(self.params.atm_max_pres), self.nlevels)
        self.pressure_profile_levels =  np.exp(press_exp)[::-1]

        # get mid point pressure between levels (i.e. get layer pressure) computing geometric
        # average between pressure at n and n+1 level
        self.pressure_profile = np.power(10, np.log10(self.pressure_profile_levels)[:-1]+
                                         np.diff(np.log10(self.pressure_profile_levels))/2.)


    def get_temperature_profile_from_file(self):

        temperature_profile = np.interp(self.pressure_profile,
                                        self.data.mixratio_pressure_file,
                                        self.data.tp_profile_file)

    def get_temperature_profile_from_param(self):

        # set initial temperature profile (mainly used for create_spectrum)
        if self.params.atm_tp_type == 'guillot':
            TP_params = [self.params.atm_tp_guillot_T_irr,
                         np.log10(self.params.atm_tp_guillot_kappa_ir),
                         np.log10(self.params.atm_tp_guillot_kappa_v1),
                         np.log10(self.params.atm_tp_guillot_kappa_v2),
                         self.params.atm_tp_guillot_alpha]
            logging.info('Set TP profile using Guillot model')
            logging.info('T_irr = %.2f' % self.params.atm_tp_guillot_T_irr)
            logging.info('kappa_ir = %.3f' % self.params.atm_tp_guillot_kappa_ir)
            logging.info('kappa_v1 = %.3f' % self.params.atm_tp_guillot_kappa_v1)
            logging.info('kappa_v2 = %.3f' % self.params.atm_tp_guillot_kappa_v2)
            logging.info('alpha = %.3f' % self.params.atm_tp_guillot_alpha)
        elif self.params.atm_tp_type == '2point':
            TP_params = [self.params.atm_tp_2point_T_surf,
                         self.params.atm_tp_2point_T_trop_diff,
                         self.params.atm_tp_2point_P_trop]
            logging.info('Set TP profile using 2point model')
            logging.info('T_surface = {0}'.format(self.params.atm_tp_2point_T_surf))
            logging.info('T_tropo_diff = {0}'.format(self.params.atm_tp_2point_T_trop_diff))
            logging.info('P_tropo = {0}'.format(self.params.atm_tp_2point_P_trop))
        elif self.params.atm_tp_type == 'Npoint':
            TP_params = [self.params.atm_tp_Npoint_T_list,
                         self.params.atm_tp_Npoint_P_list,
                         self.params.atm_tp_Npoint_smooth]
            logging.info('Set TP profile using Npoint model')
            logging.info('T_nodes = {0}'.format(self.params.atm_tp_Npoint_T_list))
            logging.info('P_nodes = {0}'.format(self.params.atm_tp_Npoint_P_list))
            logging.info('TP smoothing range = {0}'.format(self.params.atm_tp_Npoint_smooth))
        else: # if not guillot, always assume isothermal in forward model (i.e. create_spectrum)
            TP_params = [self.params.atm_tp_iso_temp]
            logging.info('Set isothermal profile with T = %.2f' % self.params.atm_tp_iso_temp)

        self.set_TP_profile(profile=self.params.atm_tp_type)

        return np.asarray(self.TP_profile(TP_params), order='C')

    # altitude, gravity and scale height profile
    def set_altitude_gravity_scaleheight_profile(self):

        # build the altitude profile from the bottom up

        H = np.zeros(self.nlayers)
        g = np.zeros(self.nlayers)
        z = np.zeros(self.nlayers)

        g[0] = (G * self.planet_mass) / (self.planet_radius**2) # surface gravity (0th layer)
        H[0] = (KBOLTZ*self.temperature_profile[0])/(self.mu_profile[0]*g[0]) # scaleheight at the surface (0th layer)

        for i in range(1, self.nlayers):
            deltaz = (-1.)*H[i-1]*np.log(self.pressure_profile_levels[i]/self.pressure_profile_levels[i-1])
            z[i] = z[i-1] + deltaz # altitude at the i-th layer

            with np.errstate(over='ignore'):
                g[i] = (G * self.planet_mass) / ((self.planet_radius + z[i])**2) # gravity at the i-th layer
            with np.errstate(divide='ignore'):
                H[i] = (KBOLTZ*self.temperature_profile[i])/(self.mu_profile[i]*g[i])

        self.altitude_profile = z
        self.scaleheight_profile = H
        self.gravity_profile = g

    # set the density profile
    def set_density_profile(self, T=None, P=None):
        self.density_profile = (self.pressure_profile)/(KBOLTZ*self.temperature_profile)

    # calculates non-linear sampling grid for TP profile from provided covariance
    def get_TP_sample_grid(self, covmat, delta=0.02):

        Pindex = [0]
        for i in range(self.nlayers-1):
            if np.abs(covmat[0,i] - covmat[0,Pindex[-1]]) >= delta:
                Pindex.append(i)
            elif (i - Pindex[-1]) >= 10:
                Pindex.append(i)
        Pindex.append(self.nlayers-1)

        self.P_index = Pindex
        self.P_sample   = self.pressure_profile[Pindex]

    # get sigma array from data.sigma_dict. Interpolate sigma array to pressure profile
    def get_sigma_array(self, nthreads=1):

        logging.info('Interpolate sigma array to pressure profile')

        pressure_profile_bar = self.pressure_profile/1e5

        sigma_array = np.zeros((self.nactivegases, len(self.pressure_profile), len(self.data.sigma_dict['t']),
                                self.int_nwngrid))

        for mol_idx, mol_val in enumerate(self.active_gases):

            sigma_in = self.data.sigma_dict['xsecarr'][mol_val]
            sigma_in_cut = sigma_in[:,:,self.int_wngrid_idxmin:self.int_wngrid_idxmax]

            if  nthreads <= 1:  # run interpolation on the current thread
                for pressure_idx, pressure_val in enumerate(pressure_profile_bar):
                    for temperature_idx, temperature_val in enumerate(self.data.sigma_dict['t']):
                        for wno_idx, wno_val in enumerate(self.int_wngrid):
                            sigma_array[mol_idx, pressure_idx, temperature_idx, wno_idx] = \
                                np.interp(pressure_val, self.data.sigma_dict['p'],
                                          sigma_in_cut[:,temperature_idx,wno_idx])

            else:

                # run interpolation with multiple threads. threads are used to compute separate pressure levels.
                # Run simultaneously up to NTHREADS pressure levels.

                logging.info('Running inteprolation with  %i threads' % nthreads)

                nruns = int(self.nlayers/nthreads)
                remainder = self.nlayers % nthreads
                for round_idx in range(nruns+1):
                    gotQueues = dict()
                    def empty_queues(jobs):
                        for i, job in enumerate(jobs):
                            if not queues[i].empty():
                                if i in gotQueues:
                                    gotQueues[i] += queues[i].get()
                                else:
                                    gotQueues[i] = queues[i].get()
                    level_min_idx = nthreads*round_idx
                    if  round_idx == nruns:
                        if remainder == 0:
                            break
                        level_max_idx = self.nlayers
                    else:
                        level_max_idx = nthreads*(round_idx+1)
                    pressure_profile_bar = self.pressure_profile[level_min_idx:level_max_idx]/1e5
                    queues = [Queue() for pressure_idx, pressure_val in enumerate(pressure_profile_bar)]
                    jobs = [MultiThread_get_sigma_array(queues[pressure_idx],
                                                         pressure_idx,
                                                         pressure_profile_bar,
                                                         self.nactivegases,
                                                         self.data.sigma_dict,
                                                         sigma_in_cut,
                                                         self.active_gases,
                                                         self.int_wngrid,
                                                         self.int_nwngrid,
                                                         self.int_wngrid_idxmin,
                                                         self.int_wngrid_idxmax) for pressure_idx, pressure_val in enumerate(pressure_profile_bar)]
                    for job in jobs:
                        job.start()
                    while any([jj.is_alive() for jj in jobs]):
                        time.sleep(0.01)  # 0.01 Wait a while before next update. Slow down updates for really long runs.
                        empty_queues(jobs)
                    for pressure_idx, pressure_val in enumerate(pressure_profile_bar):
                        sigma_array[mol_idx, round_idx*nthreads+pressure_idx, :, :] = gotQueues[pressure_idx]

        return sigma_array

    def get_ktables_array(self, nthreads=1):

        logging.info('Interpolate ktables array to pressure profile')

        pressure_profile_bar = self.pressure_profile/1e5

        kcoeff_array = np.zeros((self.nactivegases, len(self.pressure_profile),
                                 len(self.data.ktable_dict['t']), self.int_nwngrid,
                                 len(self.data.ktable_dict['weights'])))

        for mol_idx, mol_val in enumerate(self.active_gases):

            ktable_in = self.data.ktable_dict['kcoeff'][mol_val]
            if  nthreads <= 1:
                for pressure_idx, pressure_val in enumerate(pressure_profile_bar): # loop through input pressure grid
                    for temperature_idx, temperature_val in enumerate(self.data.ktable_dict['t']):# loop through temps
                        for wno_idx, wno_val in enumerate(self.int_wngrid): # loop through bins
                            for kc_idx, kc_val in enumerate(self.data.ktable_dict['weights']): # loop through k-coeff (~20)
                                # store kcoeff array
                                # pre-interpolate in pressure and slice in wavenumber to internal grid
                                kcoeff_array[mol_idx, pressure_idx, temperature_idx, wno_idx, kc_idx] = \
                                    np.interp(pressure_val, self.data.ktable_dict['p'],
                                              ktable_in[:,:,self.int_wngrid_idxmin:self.int_wngrid_idxmax,:]\
                                                       [:,temperature_idx,wno_idx,kc_idx])
            else:



                # run interpolation with multiple threads. threads are used to compute separate pressure levels.
                # Run simultaneously up to NTHREADS pressure levels.

                logging.info('Running inteprolation with  %i threads' % nthreads)

                nruns = int(self.nlayers/nthreads)
                remainder = self.nlayers % nthreads
                for round_idx in range(nruns+1):
                    gotQueues = dict()
                    def empty_queues(jobs):
                        for i, job in enumerate(jobs):
                            if not queues[i].empty():
                                if i in gotQueues:
                                    gotQueues[i] += queues[i].get()
                                else:
                                    gotQueues[i] = queues[i].get()
                    level_min_idx = nthreads*round_idx
                    if  round_idx == nruns:
                        if remainder == 0:
                            break
                        level_max_idx = self.nlayers
                    else:
                        level_max_idx = nthreads*(round_idx+1)
                    pressure_profile_bar = self.pressure_profile[level_min_idx:level_max_idx]/1e5
                    queues = [Queue() for pressure_idx, pressure_val in enumerate(pressure_profile_bar)]
                    jobs = [MultiThread_get_ktables_array(queues[pressure_idx],
                                                         pressure_idx,
                                                         pressure_profile_bar,
                                                         self.nactivegases,
                                                         self.data.ktable_dict,
                                                         ktable_in,
                                                         self.active_gases,
                                                         self.int_wngrid,
                                                         self.int_nwngrid,
                                                         self.int_wngrid_idxmin,
                                                         self.int_wngrid_idxmax) for pressure_idx, pressure_val in enumerate(pressure_profile_bar)]
                    for job in jobs:
                        job.start()
                    while any([jj.is_alive() for jj in jobs]):
                        time.sleep(0.01)  # Wait a while before next update. Slow down updates for really long runs.
                        empty_queues(jobs)
                    for pressure_idx, pressure_val in enumerate(pressure_profile_bar):
                        kcoeff_array[mol_idx, round_idx*nthreads+pressure_idx, :, :] = gotQueues[pressure_idx]

        return kcoeff_array

    def get_sigma_rayleigh_array(self):
        logging.info('Interpolate Ryleigh sigma array to pressure profile')
        sigma_rayleigh_array = np.zeros((self.nactivegases+self.ninactivegases, self.int_nwngrid))
        for mol_idx, mol_val in enumerate(self.active_gases+self.inactive_gases):
            sigma_rayleigh_array[mol_idx,:] =  self.data.sigma_rayleigh_dict[mol_val]\
                [self.int_wngrid_idxmin:self.int_wngrid_idxmax]
        return sigma_rayleigh_array
    
    def get_sigma_mie_array(self):
        ''' 
        Mie approximation, replaces rayleigh scattering. 
        Formalism taken from: Lee et al. 2013, ApJ, 778, 97
        '''
        if self.params.atm_mie_flat:
            sigma_mie = np.ones(len(self.int_wngrid))*self.mie_f
        else:
            wltmp = self.int_wlgrid #getting wavelength grid 
            a = self.mie_r
            wltmp *= 1e-4 #microns to cm 
            a *= 1e-4 #microns to cm
     
            x = 2.0 * np.pi * a/ wltmp
            Qext = 5.0 / (self.mie_q * x**(-4.0) + x**(0.2))
            sigma_mie = Qext * np.pi * (a)**(2.0) * self.mie_f 
        
        return sigma_mie 
    
    def set_sigma_mie_array(self):
        self.sigma_mie_array = self.get_sigma_mie_array()

    def get_sigma_cia_array(self):
        logging.info('Interpolate CIA sigma array to pressure profile')
        sigma_cia_array = np.zeros((len(self.params.atm_cia_pairs),
                                    len(self.data.sigma_cia_dict['t']), self.int_nwngrid))
        for pair_idx, pair_val in enumerate(self.params.atm_cia_pairs):
            sigma_cia_array[pair_idx,:,:] = self.data.sigma_cia_dict['xsecarr'][pair_val]\
                [:,self.int_wngrid_idxmin:self.int_wngrid_idxmax]
        return sigma_cia_array

    def get_cia_idx(self):
        # return the gas indexes of the molecules inside the pairs
        # used to get the mixing ratios of the individual molecules in the cpp pathintegral
        # the index refers to the full array of active_gas + inactive_gas
        cia_idx = np.zeros((len(self.params.atm_cia_pairs)*2), dtype=np.int)
        if self.params.atm_cia:
            c = 0
            for pair_idx, pair_val in enumerate(self.params.atm_cia_pairs):
                for mol_idx, mol_val in enumerate(pair_val.split('-')):
                    try:
                        cia_idx[c] = self.active_gases.index(mol_val)
                    except ValueError:
                        cia_idx[c] = self.nactivegases + self.inactive_gases.index(mol_val)
                    c += 1
        return cia_idx

    def set_ace_params(self):

        # set O, C and N abundances given metallicity (in solar units) and CO ratio
        self.O_abund_dex = np.log10(self.ace_metallicity * (10**(self.ace_O_solar-12.)))+12.
        self.N_abund_dex = np.log10(self.ace_metallicity * (10**(self.ace_N_solar-12.)))+12.
        self.C_abund_dex = self.O_abund_dex + np.log10(self.ace_co)

        # H and He don't change
        self.H_abund_dex = self.ace_H_solar
        self.He_abund_dex = self.ace_He_solar

    def set_ACE(self, mixratio_mask):

            # chemically consistent model
            vector = C.c_double*self.nlayers
            a_apt = vector()
            p_apt = vector()
            t_apt = vector()
            for i in range(self.nlayers):
               a_apt[i] = self.altitude_profile[i]/1000.
               p_apt[i] = self.pressure_profile[i]/1.e5
               t_apt[i] = self.temperature_profile[i]

            # y_out has shape (nlayers, 105). 105 is the total number of molecules computed
            y_out = ((C.c_double * 105) * self.nlayers)()

            self.ace_lib.ACE(C.byref(C.c_int(self.nlayers)),
                             C.byref(a_apt),
                             C.byref(p_apt),
                             C.byref(t_apt),
                             C.byref(C.c_double(self.He_abund_dex)),
                             C.byref(C.c_double(self.C_abund_dex)),
                             C.byref(C.c_double(self.O_abund_dex)),
                             C.byref(C.c_double(self.N_abund_dex)),
                             C.byref(y_out))
            ace_profiles = np.asarray(y_out)

            for mol_idx, mol_val in enumerate(self.params.atm_active_gases):
                self.active_mixratio_profile[mol_idx, :] = ace_profiles[:, self.data.ace_active_gases_idx[mol_idx]]

            for mol_idx, mol_val in enumerate(self.params.atm_inactive_gases):
                self.inactive_mixratio_profile[mol_idx, :] = ace_profiles[:, self.data.ace_inactive_gases_idx[mol_idx]]

            if isinstance(mixratio_mask, (np.ndarray, np.generic)):
                self.active_mixratio_profile[mixratio_mask, :] = 0


            del(y_out)
            del(a_apt)
            del(p_apt)
            del(t_apt)

            # couple mu to composition
            self.set_mu_profile()

            # update atmospheric params
            self.set_altitude_gravity_scaleheight_profile()


    ##########################################
    # Functions related to TP profile

    def set_TP_profile(self, profile=None):
        '''
        Decorator supplying TP_profile with correct TP profile function.
        Only the free parameters will be provided to the function after TP profile
        is set.
        '''
        if profile is not None: #implicit or explicit TP_type setting
            self.TP_type = profile

        if self.TP_type == 'isothermal':
            _profile = self._TP_isothermal
        elif self.TP_type == 'guillot':
            _profile = self._TP_guillot2010
        elif self.TP_type == 'Npoint':
            _profile = self._TP_Npoint
        elif self.TP_type == 'rodgers':
            _profile = self._TP_rodgers2000
        elif self.TP_type == 'hybrid':
            _profile = self._TP_hybrid
        elif self.TP_type == '2point':
            _profile = self._TP_2point
        elif self.TP_type == '3point':
            _profile = self._TP_3point
        else:
            logging.error('Invalid TP profile name')

        def tpwrap(fit_params, **kwargs):
            return self._TP_profile(_profile, fit_params, **kwargs)

        self.TP_profile = tpwrap  #setting TP_profile
        self.TP_setup = True #flag that allows TP profile internal things to be run once only.
        # Sets to False after execution

    def _TP_profile(self, TP_function, TP_params, **kwargs):
        '''
        TP profile decorator. Takes given TP profile function and fitting parameters
        To generate Temperature, Pressure, Column Density (T, P, X) grids

        fit_params depend on TP_profile selected.
        INDEX splits column densities from TP_profile parameters, INDEX = [Chi, TP]
        '''
      
        #generating TP profile given input TP_function
        T = TP_function(TP_params, **kwargs)

        self.TP_setup = False #running profile setup only on first call

        return T

    def _TP_isothermal(self, TP_params):
        '''
        TP profile for isothermal atmosphere. Follows the old implementation.
        '''

        T = np.zeros((self.nlayers))
        T[:] = TP_params

        return T

    def _TP_guillot2010(self, TP_params):
        '''
        TP profile from Guillot 2010, A&A, 520, A27 (equation 49)
        Using modified 2stream approx. from Line et al. 2012, ApJ, 749,93 (equation 19)

        Free parameters required:
            - T_irr    = planet equilibrium temperature (Line fixes this but we keep as free parameter)
            - kappa_ir = mean infra-red opacity
            - kappa_v1 = mean optical opacity one
            - kappa_v2 = mean optical opacity two
            - alpha    = ratio between kappa_v1 and kappa_v2 downwards radiation stream

        Fixed parameters:
            - T_int    = internal planetary heat flux (can be largely ignored. line puts it on 200K for hot jupiter)
            - P        = Pressure grid
            - g        = surface gravity
        '''

        # assigning fitting parameters
        T_irr = TP_params[0];
        kappa_ir = np.power(10, TP_params[1]);
        kappa_v1 = np.power(10, TP_params[2]);
        kappa_v2 = np.power(10, TP_params[3]);
        alpha = TP_params[4]

        planet_grav = (G * self.planet_mass) / (self.planet_radius**2) # surface gravity todo might change to full gravity_profile?
        gamma_1 = kappa_v1/(kappa_ir + 1e-10); gamma_2 = kappa_v2/(kappa_ir + 1e-10)
        tau = kappa_ir * self.pressure_profile / planet_grav

        T_int = 100 # todo internal heat parameter looks suspicious... needs looking at.

        def eta(gamma, tau):
            part1 = 2.0/3.0 + 2.0 / (3.0*gamma) * (1.0 + (gamma*tau/2.0 - 1.0) * np.exp(-1.0 * gamma * tau))
            part2 = 2.0 * gamma / 3.0 * (1.0 - tau**2/2.0) * spe.expn(2,(gamma*tau))
            return part1 + part2

        T4 = 3.0*T_int**4/4.0 * (2.0/3.0 + tau) + 3.0*T_irr**4/4.0 *(1.0 - alpha) * eta(gamma_1,tau) + 3.0 * T_irr**4/4.0 * alpha * eta(gamma_2,tau)
        T = T4**0.25

        return np.asarray(T)

    def _TP_rodgers2000(self, TP_params, h=None, covmatrix=None):
        '''
        Layer-by-layer temperature - pressure profile retrieval using dampening factor
        Introduced in Rodgers (2000): Inverse Methods for Atmospheric Sounding (equation 3.26)
        Featured in NEMESIS code (Irwin et al., 2008, J. Quant. Spec., 109, 1136 (equation 19)
        Used in all Barstow et al. papers.

        Free parameters required:
            - T  = one temperature per layer of Pressure (P)

        Fixed parameters:
            - P  = Pressure grid, fixed to self.pressure_profile
            - h  = correlation parameter, in scaleheights, Line et al. 2013 sets this to 7, Irwin et al sets this to 1.5
                  may be left as free and Pressure dependent parameter later.
        '''
        if h is None:
            h = self.params.atm_tp_corr_length

        # assigning parameters
        T_init = TP_params

#         covmatrix = np.zeros((self.nlayers,self.nlayers))
#         for i in range(self.nlayers):
#                 covmatrix[i,:] = np.exp(-1.0* np.abs(np.log(self.pressure_profile[i]/self.pressure_profile[:]))/h)

        if covmatrix is None: #if covariance not provided, generate
            if self.TP_setup: #run only once and save
                self.rodgers_covmat = np.zeros((self.nlayers,self.nlayers))
                for i in range(self.nlayers):
                    self.rodgers_covmat[i,:] = np.exp(-1.0* np.abs(np.log(self.pressure_profile[i]/self.pressure_profile[:]))/h)
        else:
            self.rodgers_covmat = covmatrix

        T = np.zeros((self.nlayers)) #temperature array

        #correlating temperature grid with covariance matrix
        for i in range(self.nlayers):
#             covmat  = np.exp(-1.0* np.abs(np.log(self.pressure_profile[i]/self.pressure_profile[:]))/h)
            weights = self.rodgers_covmat[i,:] / np.sum(self.rodgers_covmat[i,:])
            T[i]    = np.sum(weights*T_init)

#         plt.figure(3)
#         plt.imshow(self.rodgers_covmat,origin='lower')
#         plt.colorbar()
#         plt.show()

        #sets list of ascii parameter names. This is used in output module to compile parameters.dat
        if self.TP_setup:
            self.TP_params_ascii = []
            for i in range(self.nlayers):
                self.TP_params_ascii.append('T_{0}'.format(str(i)))
        return T

    def _TP_hybrid(self, TP_params, h=None):
        '''
        Layer-by-layer temperature pressure profile. It is a hybrid between the _TP_rodgers2000 profile and
        a second (externally calculated) covariance matrix. The external covariance can be calculated
        using a maximum likelihood retrieval of the TP profile before (or similar).
        The hybrid covariance is given by Cov_hyb = (1-alpha) * Cov_TP_rodgers200 + alpha * Cov_external

        Free parameters required:
            - alpha  = weighting parameter between covariance 1 and 2
            - T      = one temperature per layer of Pressure (P)

        Fixed parameters:
            - P  = Pressure grid, fixed to self.pressure_profile
            - h  = correlation parameter, in scaleheights, Line et al. 2013 sets this to 7, Irwin et al sets this to 1.5
                  may be left as free and Pressure dependent parameter later.
        '''
        if h is None:
            h = self.params.atm_tp_corr_length

        #if self.TP_setup:
        T = np.zeros((self.nlayers)) #temperature array

        #assigning parameters
        alpha  = TP_params[0]

        T_init = TP_params[1:]
        covmatrix = self.hybrid_covmat

        #interpolating fitting temperatures to full grid
        T_interp = np.interp(np.log(self.pressure_profile[::-1]),np.log(self.P_sample[::-1]),T_init[::-1])[::-1]

        if self.TP_setup: #run only once and save
            self.rodgers_covmat = np.zeros((self.nlayers,self.nlayers))
            for i in range(self.nlayers):
                self.rodgers_covmat[i,:] = np.exp(-1.0* np.abs(np.log(self.pressure_profile[i]/self.pressure_profile[:]))/h)

        T[:] = 0.0 #temperature array
        cov_hybrid = (1.0-alpha) * self.rodgers_covmat + alpha * covmatrix

        #correlating temperature grid with covariance matrix
        for i in range(self.nlayers):
            weights = cov_hybrid[i,:] / np.sum(cov_hybrid[i,:])
            T[i]    = np.sum(weights*T_interp)

#             plt.plot(cov_hybrid[i,:])
#         plt.ion()
#         plt.figure(2)
#         plt.clf()
# #         plt.plot(self.T,self.pressure_profile,'blue')
# #         plt.plot(T_interp,self.pressure_profile,'k')
#         plt.plot(T_init,self.P_sample,'ro')
#         plt.yscale('log')
# #         xlabel('Temperature')
# #         ylabel('Pressure (Pa)')
#         plt.gca().invert_yaxis()
#         plt.draw()
#
#         plt.ion()
#         plt.figure(3)
#         plt.clf()
#         plt.plot(self.T,self.pressure_profile,'blue')
# #         plt.plot(T_interp,self.pressure_profile,'k')
# #         plt.plot(T_init,self.P_sample,'ro')
#         plt.yscale('log')
# #         xlabel('Temperature')
# #         ylabel('Pressure (Pa)')
#         plt.gca().invert_yaxis()
#         plt.draw()


#         plt.figure(3)
#         plt.plot(weights)
#         plt.show()
#         exit()

        return T

    def set_TP_hybrid_covmat(self,covariance):
        '''
        Setting external covariance for _TP_hybrid
        '''
        self.hybrid_covmat = covariance


    def _TP_2point(self,TP_params):
        '''
        Two point TP profile. Simplest possible profile only fitting for the surface temperature and tropopause temp.
        and pressure. Atmosphere above tropopause is isothermal. Temperature is linearly interpolated in log space.
        No inversion possible here.

        Free parameters required:
            - T1 = surface temperature (at 10bar)
            - T2 = temperature at tropopause (given as difference from T1)
            - P1 = pressure at tropopause

        Fixed parameters:
            - Pressure grid (self.pressure_profile)
        '''

        maxP = np.max(self.pressure_profile)
        minP = np.min(self.pressure_profile)

        T_trop = TP_params[0] - TP_params[1]

        P_params = [maxP,TP_params[-1],minP]
        T_params = [TP_params[0],T_trop,T_trop]

        #creating linear T-P profile
        T = np.interp(np.log(self.pressure_profile[::-1]), np.log(P_params[::-1]), T_params[::-1])

        return T[::-1]


    def _TP_3point(self,TP_params):
        '''
        Same as 2point TP profile but adds one extra point between surface and troposphere

        Free parameters required:
            - T1 = surface temperature (at 10bar)
            - T2 = point 1 temperature (given as difference from T1)
            - P1 = point 1 pressure
            - T3 = point 2 temperature (given as difference from T1)
            - P2 = point 2 pressure

        Fixed parameters
            - Pressure grid (self.pressure_profile)
        '''

        maxP = np.max(self.pressure_profile)
        minP = np.min(self.pressure_profile)

        T_point1 = TP_params[0] - TP_params[1]
        T_point2 = T_point1 - TP_params[2]
        P_params = [maxP,TP_params[-2],TP_params[-1],minP]
        T_params = [TP_params[0],T_point1,T_point2, T_point2]

        #creating linear T-P profile
        T = np.interp(np.log(self.pressure_profile[::-1]), np.log(P_params[::-1]), T_params[::-1])

        return T[::-1]
    
    def _TP_Npoint(self,TP_params):
        '''
        Same as 2point TP profile but accepts N number of temperature and pressure nodes. Also includes smoothing. 
        Currently not used for fitting.
        
        Free parameters required:
            - T_array = list [] of temperature points (surface -> top)
            - P_array = list [] of corresponding pressure points (PA)
            - smooth_window = smoothing window (no. of layers)
        '''
        Tnodes = TP_params[0]
        Pnodes = TP_params[1]
        smooth_window = TP_params[2]
        
        #if first and last pressure input is -1, it replaces with P_MAX and P_MIN from 
        #atmosphere pressure grid 
        if int(Pnodes[0]) == -1 and int(Pnodes[-1]) == -1:
            Pnodes[0] = self.pressure_profile[0]
            Pnodes[-1] = self.pressure_profile[-1]          
        
        TP = np.interp((np.log(self.pressure_profile[::-1])), np.log(Pnodes[::-1]), Tnodes[::-1])
        #smoothing T-P profile
        wsize = self.nlayers*(smooth_window/100.0)
        if (wsize %2 == 0):
            wsize += 1
        TP_smooth = movingaverage(TP,wsize)
        border = np.int((len(TP) - len(TP_smooth))/2)
        
        #set atmosphere object
        foo = TP[::-1]
        foo[border:-border] = TP_smooth[::-1]

        return np.copy( foo , order='C')


# Additional classes to manage multithreading tasks

class MultiThread_get_sigma_array(Process):

    # Interpolation of sigma_array to pressure profile

    def __init__(self,
                 queue,
                 pressure_idx,
                 pressure_profile_bar,
                 nactivegases,
                 sigma_dict,
                 sigma_in_cut,
                 active_gases,
                 int_wngrid,
                 int_nwngrid,
                 int_wngrid_idxmin,
                 int_wngrid_idxmax):

        Process.__init__(self)
        self.queue = queue
        self.pressure_idx = pressure_idx
        self.pressure_profile_bar = pressure_profile_bar
        self.nactivegases = nactivegases
        self.sigma_dict = sigma_dict
        self.sigma_in_cut = sigma_in_cut
        self.active_gases = active_gases
        self.int_wngrid = int_wngrid
        self.int_nwngrid = int_nwngrid
        self.int_wngrid_idxmin = int_wngrid_idxmin
        self.int_wngrid_idxmax = int_wngrid_idxmax
        return

    def run(self):
        pressure_val = self.pressure_profile_bar[self.pressure_idx]
        if cythonised:
            sigma_array_partial = cy_fun.get_sigma_array_interp(self.sigma_dict['t'],
                                                                self.sigma_dict['p'],
                                                                pressure_val,
                                                                self.sigma_in_cut,
                                                                self.int_wngrid,
                                                                self.int_wngrid_idxmin,
                                                                self.int_wngrid_idxmax)
        else:
            sigma_array_partial = np.zeros((len(self.sigma_dict['t']), self.int_nwngrid))
            for temperature_idx, temperature_val in enumerate(self.sigma_dict['t']):
                for wno_idx, wno_val in enumerate(self.int_wngrid):
                    sigma_array_partial[temperature_idx, wno_idx] = \
                        np.interp(pressure_val, self.sigma_dict['p'], self.sigma_in_cut[:,temperature_idx,wno_idx])
        self.queue.put(sigma_array_partial)
        return

class MultiThread_get_ktables_array(Process):

    # Interpolation of sigma_array to pressure profile

    def __init__(self,
                 queue,
                 pressure_idx,
                 pressure_profile_bar,
                 nactivegases,
                 ktable_dict,
                 ktable_in,
                 active_gases,
                 int_wngrid,
                 int_nwngrid,
                 int_wngrid_idxmin,
                 int_wngrid_idxmax):

        Process.__init__(self)
        self.queue = queue
        self.pressure_idx = pressure_idx
        self.pressure_profile_bar = pressure_profile_bar
        self.nactivegases = nactivegases
        self.ktable_dict = ktable_dict
        self.ktable_in = ktable_in
        self.active_gases = active_gases
        self.int_wngrid = int_wngrid
        self.int_nwngrid = int_nwngrid
        self.int_wngrid_idxmin = int_wngrid_idxmin
        self.int_wngrid_idxmax = int_wngrid_idxmax
        return

    def run(self):
        pressure_val = self.pressure_profile_bar[self.pressure_idx]
        kcoeff_array_partial = np.zeros((len(self.ktable_dict['t']),self.int_nwngrid, len(self.ktable_dict['weights'])))
        for temperature_idx, temperature_val in enumerate(self.ktable_dict['t']):# loop through temps
            for wno_idx, wno_val in enumerate(self.int_wngrid): # loop through bins
                for kc_idx, kc_val in enumerate(self.ktable_dict['weights']): # loop through k-coeff (~20)
                    # store kcoeff array
                    # pre-interpolate in pressure and slice in wavenumber to internal grid
                    kcoeff_array_partial[temperature_idx, wno_idx, kc_idx] = \
                        np.interp(pressure_val, self.ktable_dict['p'],
                                  self.ktable_in[:,:,self.int_wngrid_idxmin:self.int_wngrid_idxmax,:]\
                                                [:,temperature_idx,wno_idx,kc_idx])

        self.queue.put(kcoeff_array_partial)
        return