'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Data class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries     
import numpy
import os
import glob
import pickle
import logging
import numpy as np

from scipy.interpolate import interp1d

from library_constants import *
from library_general import *
from library_emission import *

import license
from license import *

import matplotlib.pylab as plt

try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
except:
    MPIrank = 0
    pass

class data(object):

    def __init__(self, params):

        logging.info('Initialising data object')

        #checking taurex license
        license_manager().run()
        
        self.params = params

        # compile shared libraries
        if self.params.gen_compile_cpp:
            self.compile_shared_libs()

        # reading in atmospheric profile file from Venot
        if self.params.ven_load:
            self.load_venot_model()

        # load ace specific parameters
        if self.params.gen_ace:
            self.load_ace_params()

        # load observed input spectrum
        self.load_input_spectrum()

        # Preload the cross sections
        self.load_sigma_dict()

        # Preload rayleigh cross secitons
        self.load_sigma_rayleigh_dict()

        # Preload CIA cross sections
        self.load_sigma_cia_dict()

        # Create wavenumber grid of internal model
        self.load_wavenumber_grid()

        # reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.star_sed = self.get_star_SED()

        logging.info('Data object initialised')

    def compile_shared_libs(self):

        logging.info('Compiling shared libraries')

        if  self.params.gen_type == 'transmission':
            os.system('rm library/ctypes_pathintegral_transmission.so')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_transmission.so library/ctypes_pathintegral_transmission.cpp')
        elif  self.params.gen_type == 'emission':
            os.system('rm library/ctypes_pathintegral_emission.so')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_emission.so library/ctypes_pathintegral_emission.cpp')
        if  self.params.gen_ace:
            #os.system('rm library/ACE/ACE.so')
            os.system('gfortran -shared -fPIC  -o library/ACE/ACE.so library/ACE/Md_ACE.f90 library/ACE/Md_Constantes.f90 '
                      'library/ACE/Md_Types_Numeriques.f90 library/ACE/Md_Utilitaires.f90 library/ACE/Md_numerical_recipes.f90')

    def load_ace_params(self):

        logging.info('Loading ace specific parameters')

        # get molcules list from composes.dat
        self.ace_molecules = []
        self.ace_molecules_mu = []
        with open("library/ACE/composes.dat", "r") as textfile:
            for line in textfile:
                sl = line.split()
                self.ace_molecules.append(sl[1])
                self.ace_molecules_mu.append(float(sl[2]))

        # include He, H and N as inactive gases
        self.params.atm_inactive_gases = ['H2', 'HE', 'N2']
        self.params.atm_inactive_gases_mixratios = [0.85, 0.15, 0]
        self.params.atm_inactive_gases_idx = [0, 0, 0] # will be set below

        # include only molecules for which we have cross sections
        self.params.atm_active_gases = []
        self.params.atm_active_gases_mixratios = []
        self.ace_active_gases_idx = []
        self.ace_inactive_gases_idx = [0,0,0]

        for mol_idx, mol_val in enumerate(self.ace_molecules):
            molpath = os.path.join(self.params.in_xsec_path, '%s.db' % mol_val)
            if os.path.isfile(molpath):
                self.params.atm_active_gases.append(mol_val)
                self.params.atm_active_gases_mixratios.append(0)
                self.ace_active_gases_idx.append(mol_idx)

            if mol_val.upper() == 'H2':
                self.ace_inactive_gases_idx[0] = mol_idx
            if mol_val.upper() == 'HE':
                self.ace_inactive_gases_idx[1] = mol_idx
            if mol_val.upper() == 'N2':
                self.ace_inactive_gases_idx[2] = mol_idx


        logging.info('ace active absorbers: %s' % self.params.atm_active_gases)
        logging.info('ace active absorbers idx: %s' % self.ace_active_gases_idx)
        logging.info('ace inctive absorbers: %s' % self.params.atm_inactive_gases)
        logging.info('ace inactive absorbers idx: %s' % self.ace_inactive_gases_idx)

    def load_input_spectrum(self):

        # set observed spectrum specific variables (only if spectrum is provided)
        if self.params.in_spectrum_file != 'False':

            # read spectrum from file
            logging.info('Reading spectrum from file: %s ' % self.params.in_spectrum_file)
            spectrum = np.loadtxt(self.params.in_spectrum_file)
            spectrum = spectrum[spectrum[:,0].argsort(axis=0)]
            self.obs_spectrum = spectrum
            self.obs_wlgrid = self.obs_spectrum[:,0] # grid in micron
            self.obs_wngrid = 10000./self.obs_spectrum[:,0] # grid in wavenumbers
            self.obs_nwlgrid = len(self.obs_spectrum[:,0]) # number of datapoints in spectrum
            self.obs_binwidths = self.obs_spectrum[:,3]   if shape(self.obs_spectrum)[1] == 4 else None # bin widths
        else:

            # no input spectrum provided. If running only forward model, input spectrum can be omitted
            logging.info('No input spectrum provided')
            self.obs_spectrum = False

    def load_venot_model(self):

            logging.info('Load atmospheric profile from external files. Format: Venot photochemical models')

            # altitude (km), pressure (mbar), temperature (K) profile
            apt = np.loadtxt(self.params.ven_TP_profile_path)
            self.ven_altitude = apt[:,0]/1000. # convert km to m
            self.ven_pressure = apt[:,1]*100. # convert mbar to pascal
            self.ven_temperature = apt[:,2]
            self.ven_altitude_int = interp1d(self.ven_pressure, self.ven_altitude)
            self.ven_temperature_int = interp1d(self.ven_pressure, self.ven_temperature)

            # mixing ratio profiles

            # extract list of molecules from 'fractions_molaires' file
            with open(self.params.ven_mol_profile_path, 'r') as f:
                self.ven_molecules = [x.upper() for x in f.readline().split()]
                self.ven_molweight = [float(x) for x in f.readline().split()] # molecular weight in AMU
                table = np.loadtxt(self.params.ven_mol_profile_path, skiprows=2)
                self.ven_molprof_altitude = table[:,0]/1000. # convert km to m
                self.ven_molprof_pressure = table[:,1]*100. # convert mbar to pascal
                self.ven_molprof_mixratios = table[:,2:] # mixing ratios referring to self.ven_molecules
                self.ven_molprof_altitude_int = interp1d(self.ven_molprof_pressure, self.ven_molprof_altitude)
                self.ven_molprof_mixratios_int = [interp1d(self.ven_molprof_pressure, self.ven_molprof_mixratios[:,i]) for i in xrange(np.shape(self.ven_molprof_mixratios)[1])]
                logging.info('Atmospheric pressure boundaries from chemical model: %f-%f' % (np.min(self.ven_molprof_pressure), np.max(self.ven_molprof_pressure)))

            # determine list of molecules with cross sections. exclude the others.
            self.ven_active_gases = []
            self.ven_active_gases_idx = []
            self.ven_noxsec = []

            for mol_idx, mol_val in enumerate(self.ven_molecules):
                molpath = os.path.join(self.params.in_xsec_path, '%s.db' % mol_val)
                if os.path.isfile(molpath) and not mol_val in self.params.ven_exclude_mol:
                    self.ven_active_gases.append(mol_val)
                    self.ven_active_gases_idx.append(mol_idx)
                else:
                    self.ven_noxsec.append(mol_val)

            logging.info('Molecules with available cross sections: %s' %  self.ven_active_gases)
            logging.info('Molecules with available cross sections idx: %s' %  self.ven_active_gases_idx)
            logging.info('Excluded molecules: %s' %  self.ven_noxsec)

            # determine list of inactive gases
            gases = ['H2', 'HE', 'N2']
            self.ven_inactive_gases = []
            self.ven_inactive_gases_idx = []
            for gasname in gases:
                if gasname in self.ven_molecules:
                    self.ven_inactive_gases.append(gasname)
                    self.ven_inactive_gases_idx.append(self.ven_molecules.index(gasname))
            logging.info('Inactive gases: %s' %  self.ven_inactive_gases)

    def load_sigma_dict(self):

        # preload the absorption cross sections for the input molecules
        if self.params.ven_load:
            molecules = self.ven_active_gases
        else:
            molecules = self.params.atm_active_gases

        sigma_dict = {}
        sigma_dict['xsecarr'] = {}

        for mol_idx, mol_val in enumerate(molecules):

            molpath = os.path.join(self.params.in_xsec_path, '%s.db' % mol_val)
            if not os.path.isfile(molpath):
                logging.error('There is no cross section for %s. Path: %s ' % (mol_val, molpath))
                exit()
            else:
                sigma_tmp = pickle.load(open(molpath)) # load pickled cross section array for molecule mol_val

                # check that the wavenumber, temperature and pressure grid are the same for all cross sections
                if mol_idx > 0:
                    if np.unique(sigma_tmp['wno'] - wno) != 0:
                        logging.error('The cross section of %s has a different wavenumber grid than the one of %s' % (molecules[mol_idx], molecules[mol_idx-1]))
                        exit()
                    if np.unique(sigma_tmp['t'] - t) != 0:
                        logging.error('The cross section of %s has a different temperature grid than the one of %s' % (molecules[mol_idx], molecules[mol_idx-1]))
                        exit()
                    if np.unique(sigma_tmp['p'] - p) != 0:
                        logging.error('The cross section of %s has a different pressure grid than the one of %s' % (molecules[mol_idx], molecules[mol_idx-1]))
                        exit()

                t = sigma_tmp['t']
                p = sigma_tmp['p']
                wno = sigma_tmp['wno']

                if mol_idx == 0:
                #     if np.min(wno) > np.min(self.int_wngrid) or np.max(wno) < np.max(self.int_wngrid):
                #         logging.error('Internal wavenumber grid overflow')
                #         logging.error('Internal wavenumber grid range: %f - %f' % (np.min(self.int_wngrid), np.max(self.int_wngrid)))
                #         logging.error('Cross section wavenumber grid range: %f - %f' % (np.min(wno), np.max(wno)))
                #         exit()
                #     if np.unique(np.diff(wno)) != self.params.in_xsec_dnu:
                #         logging.error('The resolution of the internal model (%f wno) is different from '
                #                       'the resolution of the cross sections (%f wno)' % (self.params.in_xsec_dnu, np.unique(np.diff(wno))))
                #         exit()

                    # restrict wavenumber range
                    # wno_min_idx = np.where(np.abs(wno - np.min(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]
                    # wno_max_idx = np.where(np.abs(wno - np.max(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]+1

                    # restrict temperature range

                    if self.params.in_xsec_alltemp:
                        sigma_dict['t'] = t
                        Tmax_idx = len(t)
                        Tmin_idx = 0
                    else:
                        Tmin = 0
                        Tmax = 5000
                        if self.params.ven_load:
                            Tmax = np.max(self.ven_temperature)
                            Tmin = np.min(self.ven_temperature)
                        elif self.params.downhill_run or self.params.mcmc_run or self.params.nest_run:
                            if self.params.atm_tp_type == 'isothermal':
                                Tmax = self.params.fit_tp_iso_bounds[1]
                                Tmin = self.params.fit_tp_iso_bounds[0]
                        elif self.params.atm_tp_type == 'isothermal':
                            Tmin = Tmax = self.params.atm_tp_iso_temp

                        if Tmax > np.max(t) or Tmin < np.min(t):
                            logging.warning('The atmospheric temperature profile falls outside the temperature range of the cross sections')
                            logging.warning('Internal temperature range: %i - %i' % (Tmin, Tmax))
                            logging.warning('Cross section temperature range: %i - %i' % (np.min(t), np.max(t)))

                        # get idx of the temperature boundaries in the cross sections
                        # always get the max T availble in the xsec that is lower than Tmin and the min T that is higher than Tmax

                        if Tmax > np.max(t):
                            Tmax_idx = len(t)
                        else:
                            Tmax_diff = Tmax - t # t is sigma_tmp['t']
                            Tmax_idx = len(Tmax_diff) - np.argmax(Tmax_diff[::-1][Tmax_diff<=0])

                        if Tmin < np.min(t):
                            Tmin_idx = 0
                        else:
                            Tmin_diff = Tmin - t # t is sigma_tmp['t'], the list of temperature
                            Tmin_idx = np.argmin(Tmin_diff[Tmin_diff>=0])

                        if Tmin > np.max(t) and Tmax > np.max(t):
                            sigma_dict['t'] = np.asarray([t[-1]])
                            Tmax_idx = len(t)
                        else:
                            sigma_dict['t'] = t[Tmin_idx:Tmax_idx]

                    sigma_dict['p'] = p

                    #sigma_dict['wno'] = sigma_tmp['wno'][wno_min_idx:wno_max_idx] # check: this should be identical to internal model!
                    sigma_dict['wno'] = sigma_tmp['wno'] # check: this should be identical to internal model!

                    # if np.unique(sigma_dict['wno'] - self.int_wngrid) != 0:
                    #     logging.error('Wavenumber grid of xsec is different from wn grid of internal model')
                    #     exit()


                #sigma_dict['xsecarr'][mol_val] = sigma_tmp['xsecarr'][:,Tmin_idx:Tmax_idx,wno_min_idx:wno_max_idx] / 10000. # also convert from cm^-2 to m^-2
                sigma_dict['xsecarr'][mol_val] = sigma_tmp['xsecarr'][:,Tmin_idx:Tmax_idx] / 10000. # also convert from cm^-2 to m^-2
        logging.info('Temperature range: %s' % sigma_dict['t'] )
        logging.info('Pressure range: %s ' % sigma_dict['p'])

        logging.info('The full wavenumber range is %.2f - %.2f with %i points' % (np.min(wno), np.max(wno), len(wno)))

        self.int_wngrid_full = wno # full wavenumber range
        self.int_nwngrid_full = len(wno)

        self.int_wlgrid_full = 10000./wno
        self.int_nwlgrid_full = len(wno)

        self.sigma_dict = sigma_dict

        del sigma_tmp, t, p, wno

    def load_sigma_rayleigh_dict(self):

        # precalculate rayleigh scattering cross sections

        logging.info('Compute Rayleigh scattering cross sections')

        sigma_rayleigh_dict = {}

        if self.params.ven_load:
            molecules = self.ven_active_gases + self.ven_inactive_gases
        else:
            molecules = self.params.atm_active_gases + self.params.atm_inactive_gases

        for gasname in molecules:

            gasname = gasname.upper()

            # get the refractive index. Formulae taken from Allen Astrophysical Quantities if not otherwise specified
            n_formula = True # assume we have a formula for the refractive index of gasname
            king = 1 # King correction factor
            ns = 0   # refractive index
            wn = self.int_wngrid_full # wavenumber in cm^-1

            if gasname == 'HE':
                ns = 1 + 0.01470091/(423.98-(10000./wn)**-2) # C. R. Mansfield and E. R. Peck. Dispersion of helium, J. Opt. Soc. Am. 59, 199-203 (1969)
                # king is one for He
            elif gasname == 'H2':
                ns = 1 + 13.58e-5 * (1. + 7.52e-3 / (10000./wn)**2)
                delta = 0.035 # from Morgan old code..
                king = (6.+3.*delta)/(6.-7.*delta) # Bates (1984)
            elif gasname == 'N2':
                ns = 1. + (6498.2 + (307.43305e12)/(14.4e9 - wn**2))*1.e-8 # Peck and Khanna
                king = 1.034+3.17e-12*wn**2 # Bates
            elif gasname == 'O2':
                ns = 1 + 1.181494e-4 + 9.708931e-3/(75.4-(10000./wn)**-2) #  J. Zhang, Z. H. Lu, and L. J. Wang. Appl. Opt. 47, 3143-3151 (2008)
                king = 1.096
            elif gasname == 'CO2':
                #A. Bideau-Mehu, Y. Guern, R. Abjean and A. Johannin-Gilles. Interferometric determination of the refractive index of carbon dioxide in the ultraviolet region, Opt. Commun. 9, 432-434 (1973)
                ns = 1 + 6.991e-2/(166.175-(10000./wn)**-2)+1.44720e-3/(79.609-(10000./wn)**-2)+6.42941e-5/(56.3064-(10000./wn)**-2)+5.21306e-5/(46.0196-(10000./wn)**-2)+1.46847e-6/(0.0584738-(10000./wn)**-2)
                king = 1.1364 #  Sneep & Ubachs 2005
            elif gasname == 'CH4':
                ns = 1 + 1.e-8*(46662. + 4.02*1.e-6*(1/((10000./wn)*1.e-4))**2)
            elif gasname == 'CO':
                ns = 1 + 32.7e-5 * (1. + 8.1e-3 / (10000./wn)**2)
                king = 1.016  #  Sneep & Ubachs 2005
            elif gasname == 'NH3':
                ns = 1 + 37.0e-5 * (1. + 12.0e-3 / (10000./wn)**2)
            elif gasname == 'H2O':
                ns_air = (1 + (0.05792105/(238.0185 - (10000./wn)**-2) + 0.00167917/(57.362-(10000./wn)**-2))) # P. E. Ciddor. Appl. Optics 35, 1566-1573 (1996)
                ns = 0.85 * (ns_air - 1.) + 1  # ns = 0.85 r(air) (Edlen 1966)
                delta = 0.17 # Marshall & Smith 1990
                king = (6.+3.*delta)/(6.-7.*delta)
            else:
                # this sets sigma_R to zero for all other gases
                n_formula = False
                logging.warning('There is no formula for the refractive index of %s. Cannot compute the cross section' % gasname)


            if n_formula: # only if the refractive index was computed
                Ns = 2.6867805e25 # in m^-3
                sigma =  ( 24.*np.pi**3/ (Ns**2) ) * (((ns**2 - 1.)/(ns**2 + 2.))**2) * king / ((10000./wn) * 1.e-6)**4

                # override H2 and He sigma with formulae from M Line
                if gasname == 'H2':
                    wave = (1/wn)*1E8
                    sigma = ((8.14E-13)*(wave**(-4.))*(1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4
                if gasname == 'HE':
                    wave = (1/wn)*1E8
                    sigma =  ((5.484E-14)*(wave**(-4.))*(1+(2.44E5)*(wave**(-2.))))*1E-4
                    sigma[:] = 0

                logging.info('Rayleigh scattering cross section of %s correctly computed' % (gasname))
                sigma_rayleigh_dict[gasname] = sigma
            else:
                sigma_rayleigh_dict[gasname] = np.zeros((len(wn)))

        self.sigma_rayleigh_dict = sigma_rayleigh_dict

    def load_sigma_cia_dict(self):

        # preload collision induced absorption cross sections

        sigma_dict = {}
        sigma_dict['xsecarr'] = {}

        for pair_val in self.params.atm_cia_pairs:

            cia_path = os.path.join(self.params.in_cia_path, '%s.db' % pair_val.upper())

            sigma_tmp = pickle.load(open(cia_path)) # load pickled cross section array

            t = sigma_tmp['t']
            wno = sigma_tmp['wno']

            # check cia wavenumber boundaries
            if np.min(wno) > np.min(self.int_wngrid_full) or np.max(wno) < np.max(self.int_wngrid_full):
                logging.warning('Internal wavenumber grid overflow for CIA xsec for %s' % pair_val)
                logging.warning('Internal wavenumber grid range: %f - %f' % (np.min(self.int_wngrid_full), np.max(self.int_wngrid_full)))
                logging.warning('CIA cross section wavenumber grid range: %f - %f' % (np.min(wno), np.max(wno)))
                logging.warning('Assume xsec to be zero outside the cia range')


            # restrict temperature range
            if self.params.in_xsec_alltemp:
                sigma_dict['t'] = t
                Tmax_idx = len(t)
                Tmin_idx = 0
            else:
                Tmin = 0
                Tmax = 5000
                if self.params.ven_load:
                    Tmax = np.max(self.ven_temperature)
                    Tmin = np.min(self.ven_temperature)
                elif self.params.downhill_run or self.params.mcmc_run or self.params.nest_run:
                    if self.params.atm_tp_type == 'isothermal':
                        Tmax = self.params.fit_tp_iso_bounds[1]
                        Tmin = self.params.fit_tp_iso_bounds[0]
                elif self.params.atm_tp_type == 'isothermal':
                    Tmin = Tmax = self.params.atm_tp_iso_temp

                if Tmax > np.max(t) or Tmin < np.min(t):
                    logging.warning('The atmospheric temperature profile falls outside the temperature range of the cross sections')
                    logging.warning('Internal temperature range: %i - %i' % (Tmin, Tmax))
                    logging.warning('Cross section temperature range: %i - %i' % (np.min(t), np.max(t)))

                # get idx of the temperature boundaries in the cross sections
                # always get the max T availble in the xsec that is lower than Tmin and the min T that is higher than Tmax

                if Tmax > np.max(t):
                    Tmax_idx = len(t)
                else:
                    Tmax_diff = Tmax - t # t is sigma_tmp['t']
                    Tmax_idx = len(Tmax_diff) - np.argmax(Tmax_diff[::-1][Tmax_diff<=0])

                if Tmin < np.min(t):
                    Tmin_idx = 0
                else:
                    Tmin_diff = Tmin - t # t is sigma_tmp['t'], the list of temperature
                    Tmin_idx = np.argmin(Tmin_diff[Tmin_diff>=0])

                if Tmin > np.max(t) and Tmax > np.max(t):
                    sigma_dict['t'] = np.asarray([t[-1]])
                    Tmax_idx = len(t)
                else:
                    sigma_dict['t'] = t[Tmin_idx:Tmax_idx]

            sigma_dict['wno'] = self.int_wngrid_full
            sigma_dict['xsecarr'][pair_val] = np.zeros((len(sigma_dict['t']), self.int_nwngrid_full))

            # reinterpolate cia xsec to molecule cross section grid and save
            for t_idx, t_val in enumerate(sigma_dict['t']):
                sigma_dict['xsecarr'][pair_val][t_idx,:] = np.interp(self.int_wngrid_full, wno, sigma_tmp['xsecarr'][Tmin_idx+t_idx])

            # load the sigma array in memory
            logging.info('Preload cia cross section for %s' % pair_val)

        del sigma_tmp, t, wno

        self.sigma_cia_dict = sigma_dict

    def load_wavenumber_grid(self):

        logging.info('Create wavenumber grid of internal model')

        # wavenumber grid limits of internal model
        if self.params.gen_manual_waverange or not isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # limits defined by a manual wavelength range in micron in param file
            lambdamax = self.params.gen_wavemax
            lambdamin = self.params.gen_wavemin
            #self.int_wngrid_obs_idxmin = 0
            #self.int_wngrid_obs_idxmax = len(self.int_wngrid_full)
        else:
            # limits defined by the input spectrum in micron.
            lambdamin = self.obs_wlgrid[0]
            lambdamax = self.obs_wlgrid[-1]
            # Expand to half a bin up, and half a bin down to properly model edges
            if self.obs_binwidths == None:
                # if bin widths are *not* provided in the input spectrum
                bin_up =  (self.obs_wlgrid[-1]-self.obs_wlgrid[-2])/2.
                bin_low = (self.obs_wlgrid[1]-self.obs_wlgrid[0])/2.
            else:
                # if bin widths are provided in input spectrum
                bin_up = self.obs_binwidths[-1]/2.
                bin_low = self.obs_binwidths[0]/2.
            lambdamin = self.obs_wlgrid[0] - bin_low
            lambdamax = self.obs_wlgrid[-1] + bin_up

        # convert to wavenumbers
        numin = 10000./lambdamax
        numax = 10000./lambdamin

        # find numin / numax closest to the cross section wavenumber grid (approximate numin for defect, and numax for excess)
        idx_min = np.argmin(np.abs(self.int_wngrid_full-numin))
        if numin - self.int_wngrid_full[idx_min] < 0:
            idx_min -= 1

        idx_max = np.argmin(np.abs(self.int_wngrid_full-numax))
        if numax - self.int_wngrid_full[idx_max] > 0:
            idx_max += 1

        self.int_wngrid_obs_idxmin = idx_min
        self.int_wngrid_obs_idxmax = idx_max

        self.int_wngrid_obs = self.int_wngrid_full[self.int_wngrid_obs_idxmin:self.int_wngrid_obs_idxmax]
        self.int_nwngrid_obs = len(self.int_wngrid_obs)

        logging.info('Internal wavenumber grid is %.2f - %.2f with %i points' %
                     (self.int_wngrid_full[self.int_wngrid_obs_idxmin],
                      self.int_wngrid_full[self.int_wngrid_obs_idxmax-1],
                      self.int_nwngrid_obs))

        # convert wavenumber grid to wavelenght grid
        self.int_wlgrid_obs = 10000./self.int_wngrid_obs
        self.int_nwlgrid_obs = len(self.int_wlgrid_obs)

        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # calculate spectral binning grid in wavelength space
            self.intsp_bingrid, self.intsp_bingrididx = get_specbingrid(self.obs_wlgrid, self.int_wlgrid_obs, self.obs_binwidths)
            self.intsp_nbingrid = len(self.obs_wlgrid)

            self.intsp_bingrid_full, self.intsp_bingrididx_full = get_specbingrid(self.obs_wlgrid, self.int_wlgrid_full, self.obs_binwidths)
            self.intsp_nbingrid_full = len(self.obs_wlgrid)


    def get_star_SED(self):

        # reading in phoenix spectra from folder specified in parameter file

        index = loadtxt(os.path.join(self.params.in_star_path, "SPTyp_KH.dat"), dtype='string')
        tmpind = []
        for i in range(len(index)):
            tmpind.append(float(index[i][1]))
        tmpind = sort(tmpind)
        
        # reading in stellar file index
        fileindex = glob.glob(os.path.join(self.params.in_star_path, "*.fmt"))
        
        if self.params.star_temp > max(tmpind) or self.params.star_temp < min(tmpind):
            if self.params.verbose:
                logging.warning('Stellar temp. in .par file exceeds range %.1f - %.1f K. '
                                'Using black-body approximation instead' % (min(tmpind), max(tmpind)))
            self.star_blackbody = True
            SED = black_body(self.int_wngrid_obs,self.params.star_temp) #@todo bug here? not multiplied by size of star 4piRs^2
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = find_nearest(tmpind, self.params.star_temp)
            self.star_blackbody = False
            
            for file in fileindex: #this search is explicit due to compatibility issues with Mac and Linux sorting
                if np.int(file.split('/')[-1][3:8]) == np.int(tmpselect):
                    self.SED_filename = file

            #reading in correct file and interpolating it onto self.int_wngrid_obs
            SED_raw = np.loadtxt(self.SED_filename, dtype='float', comments='#')
            SED_raw[:,1] *= 10.0  #converting from ergs to SI @todo move converting somewhere more sane 

            SED = np.interp(self.int_wlgrid_obs, SED_raw[:,0], SED_raw[:,1])

        return SED
