'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Data class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries     
import numpy
import os
import glob

try:
    import cPickle as pickle
except:
    import pickle

import logging
import numpy as np
import heapq
from scipy.interpolate import interp1d

from library_constants import *
from library_general import *
from library_emission import *

#import license
#from license import *

from matplotlib.pylab import *


class data(object):

    def __init__(self, params):

        logging.info('Initialising data object')

        # DISABLED.
        #checking taurex license
        #license_manager().run()
        
        self.params = params

        # compile shared libraries
        if self.params.gen_compile_cpp:
            self.compile_shared_libs()

        # set opacity method (ktab, xsec_sampled, xsec_highres)
        self.get_opacity_method()

        # load ace specific parameters
        if self.params.gen_ace:
            self.load_ace_params()


        # load observed input spectrum
        if self.params.mode == 'retrieval':
            self.load_input_spectrum()



        # Load external active gases profiles
        if self.params.atm_active_gases[0] == 'FILE':
            self.load_active_gases_file()

        # Load external temperature profile
        if self.params.atm_tp_type.upper() == 'FILE':
            self.load_tp_profile_file()

        # Preload the cross sections
        if self.opacity_method == 'xsec_sampled':
            self.load_sigma_sampled_dict()
        elif self.opacity_method == 'xsec_highres':
            self.load_sigma_highres_dict()
        elif self.opacity_method == 'ktables':
            self.load_ktables_dict()

        # Load wavenumber grid of internal model
        self.load_wavenumber_grid()

        # Load rayleigh cross secitons
        self.load_sigma_rayleigh_dict()

        # Load CIA cross sections
        self.load_sigma_cia_dict()
        
        # Load MIE coefficients 
        self.load_mie_indices()

        # Load Phoenix stellar model library or star blackbody
        if self.params.gen_type.upper() == 'EMISSION':
            self.load_star_SED()

        logging.info('Data object initialised')

    def compile_shared_libs(self):

        # todo check. Do we want it here?

        logging.info('Compiling shared libraries')

        if  self.params.gen_type == 'transmission':
            os.system('rm library/ctypes_pathintegral_transmission_xsec.so')
            os.system('rm library/ctyped_pathintegral_transmission_ktab.so')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_transmission_xsec.so '
                      'library/ctypes_pathintegral_transmission_xsec.cpp')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_transmission_ktab.so '
                      'library/ctypes_pathintegral_transmission_ktab.cpp')
        elif  self.params.gen_type == 'emission':
            os.system('rm library/ctypes_pathintegral_emission.so')
            os.system('rm library/ctyped_pathintegral_emission_ktab.so')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_emission.so '
                      'library/ctypes_pathintegral_emission.cpp')
            os.system('g++ -fPIC -shared -o library/ctypes_pathintegral_emission_ktab.so '
                      'library/ctypes_pathintegral_emission_ktab.cpp')
        if  self.params.gen_ace:
            #os.system('rm library/ACE/ACE.so')
            os.system('gfortran -shared -fPIC  -o library/ACE/ACE.so library/ACE/Md_ACE.f90 '
                      'library/ACE/Md_Constantes.f90 '
                      'library/ACE/Md_Types_Numeriques.f90 library/ACE/Md_Utilitaires.f90 '
                      'library/ACE/Md_numerical_recipes.f90')
            
        if self.params.atm_mie and self.params.atm_mie_type == 'bh':
            os.system('gcc -fPIC -shared library/MIE/bhmie_lib.c library/MIE/complex.c library/MIE/nrutil.c -o library/MIE/bhmie_lib.so')

    def get_opacity_method(self):

        if self.params.in_opacity_method in ['xsec_sampled', 'xsec_lowres', 'xsec']:
            self.opacity_method = 'xsec_sampled'
        elif self.params.in_opacity_method == 'xsec_highres':
            self.opacity_method = 'xsec_highres'
        elif self.params.in_opacity_method in ['ktab', 'ktable', 'ktables']:
            self.opacity_method = 'ktables'
        else:
            logging.error('You need to select an opacity calculation method. See parameter General->opacity_method')
            exit()

    def load_ace_params(self):

        logging.info('Loading ace specific parameters')

        # get molcules list from composes.dat
        self.ace_molecules = []
        self.ace_molecules_mu = []
        with open('library/ACE/composes.dat', 'r') as textfile:
            for line in textfile:
                sl = line.split()
                self.ace_molecules.append(sl[1])
                self.ace_molecules_mu.append(float(sl[2]))

        # include only molecules for which we have cross sections
        self.params.atm_active_gases = []
        self.params.atm_active_gases_mixratios = []
        self.ace_active_gases_idx = []
        self.ace_inactive_gases_idx = [0, 0, 0]

        for mol_idx, mol_val in enumerate(self.ace_molecules):
            if self.opacity_method in ['xsec_sampled', 'xsec_highres']:
                opacity_path = self.params.in_xsec_path
            elif self.opacity_method == 'ktables':
                opacity_path = self.params.in_ktab_path

            molpath = insensitive_glob(os.path.join(opacity_path, '%s_*' % mol_val))+\
                      insensitive_glob(os.path.join(opacity_path, '%s.*' % mol_val))
            if len(molpath)>0:
                molpath = molpath[0]
                if os.path.isfile(molpath):
                    self.params.atm_active_gases.append(mol_val)
                    self.params.atm_active_gases_mixratios.append(0)
                    self.ace_active_gases_idx.append(mol_idx)

            # cannot find molecule, maybe it's an inactive gas?
            elif mol_val.upper() == 'H2':
                    self.ace_inactive_gases_idx[0] = mol_idx
            elif mol_val.upper() == 'HE':
                    self.ace_inactive_gases_idx[1] = mol_idx
            elif mol_val.upper() == 'N2':
                    self.ace_inactive_gases_idx[2] = mol_idx
            else:
                continue

        logging.info('ace active absorbers: %s' % self.params.atm_active_gases)
        logging.info('ace active absorbers idx: %s' % self.ace_active_gases_idx)
        logging.info('ace inctive absorbers: %s' % self.params.atm_inactive_gases)
        logging.info('ace inactive absorbers idx: %s' % self.ace_inactive_gases_idx)

        if len(self.params.atm_active_gases) == 0:
            logging.error('There are no enough active absorbers. Check that you have the right cross sections or ktables')
            exit()

    def load_input_spectrum(self):

        # set observed spectrum specific variables (only if spectrum is provided)
        if self.params.in_spectrum_file in ['False', 'None', '']:
            logging.error('You need to specify an input spectrum when running a retrieval. See Input->spectrum_file')
            exit()
        elif not os.path.isfile(self.params.in_spectrum_file):
            logging.error('The input spectrum file specified in Input->spectrum_file does not exist. Path: %s'
                          % self.params.in_spectrum_file)
            exit()
        else:

            # read spectrum from file
            logging.info('Reading spectrum from file: %s ' % self.params.in_spectrum_file)

            # input spectrum must be wavelength (micron)
            spectrum = np.loadtxt(self.params.in_spectrum_file)
            spectrum = spectrum[spectrum[:,0].argsort(axis=0)[::-1]] # sort in wavenumber

            self.obs_spectrum = spectrum
            self.obs_wlgrid = self.obs_spectrum[:,0]
            self.obs_wngrid = 10000./self.obs_wlgrid
            self.obs_nwlgrid = len(self.obs_spectrum[:,0]) # number of datapoints in spectrum

            # spectral bins widths and edges
            if shape(self.obs_spectrum)[1] == 4:

                # bins given as input
                self.obs_binwidths = self.obs_spectrum[:,3]

                # calculate bin edges
                bin_edges = []
                for i in range(len(self.obs_wlgrid)):
                    bin_edges.append(self.obs_wlgrid[i]-self.obs_binwidths[i]/2.)
                    bin_edges.append(self.obs_wlgrid[i]+self.obs_binwidths[i]/2.)
                self.obs_binedges = bin_edges

            else:

                # bins widths and edges need to be extrapolated
                bin_edges =[]
                wavegrid = self.obs_wlgrid
                bin_edges.append(self.obs_wlgrid[0]- (self.obs_wlgrid[1]-self.obs_wlgrid[0])/2.0) #first bin edge
                for i in range(len(self.obs_wlgrid)-1):
                    bin_edges.append(self.obs_wlgrid[i]+(self.obs_wlgrid[i+1]-self.obs_wlgrid[i])/2.0)
                bin_edges.append((self.obs_wlgrid[-1]-self.obs_wlgrid[-2])/2.0 + self.obs_wlgrid[-1]) #last bin edge
                bin_widths = np.abs(np.diff(bin_edges))

                self.obs_binwidths = bin_widths
                self.obs_binedges = bin_edges

            logging.info('Input spactrum wavenumber range (including bin widths): %.3f, %.3f' % (self.obs_binedges[1],
                                                                                                 self.obs_binedges[-1]))


    def load_ktables_dict(self):

        # preload the k tables for the input molecules

        if self.params.atm_active_gases[0] == 'FILE':
            molecules = self.active_gases_file
        else:
            molecules = self.params.atm_active_gases

        ktable_dict = {}
        ktable_dict['kcoeff'] = {}

        # loop through molecules
        for mol_idx, mol_val in enumerate(molecules):

            # check that ktable for given molecule exists
            molpath = insensitive_glob(os.path.join(self.params.in_ktab_path, '%s_*' % mol_val))+ \
                      insensitive_glob(os.path.join(self.params.in_ktab_path, '%s.*' % mol_val))
            if len(molpath) > 0:
                molpath = molpath[0]
            else:
                logging.error('There is no ktable for %s. Path: %s ' % (mol_val, molpath))
                exit()

            # load ktable
            try:
                ktable = pickle.load(open(molpath, 'rb'), encoding='latin1') # python 3
            except:
                ktable = pickle.load(open(molpath)) # python 2

            if mol_idx > 0:
                # check that ktables are all consistent with each others
                if not np.array_equal(t, ktable['t'])  or not np.array_equal(p, ktable['p'])  or \
                    not np.array_equal(samples, ktable['samples']) or not np.array_equal(weights, ktable['weights']) or\
                    not np.array_equal(bin_centers, ktable['bin_centers']) or \
                    not np.array_equal(bin_edges, ktable['bin_edges']) or \
                    ngauss != ktable['ngauss'] or method != ktable['method'] or wnrange != ktable['wnrange'] or \
                    wlrange != ktable['wlrange']  or wnrange != ktable['wnrange']:  #or \
                    #resolution != ktable['resolution']:
                    logging.error('Something is wrong with this ktable: %s. Check that the format of all '
                                  'ktables is consistent' % molpath)
                    exit()
            t = ktable['t']
            p = ktable['p']
            method = ktable['method']
            ngauss = ktable['ngauss']
            samples = ktable['samples']
            weights = ktable['weights']
            bin_centers = ktable['bin_centers']
            bin_edges = ktable['bin_edges']
            wnrange = ktable['wnrange']
            wlrange = ktable['wlrange']
            #resolution = ktable['resolution']
            kcoeff = ktable['kcoeff']


            if mol_idx == 0:

                # first molecule, set temperature and pressure lists in ktable_dict

                # restrict temperature range
                T_list, Tmin_idx, Tmax_idx = self.get_temp_range_idx(t)
                ktable_dict['t'] = T_list.astype(float)
                ktable_dict['p'] = p.astype(float)
                ktable_dict['bin_centers'] = bin_centers
                ktable_dict['bin_edges'] = bin_edges
                ktable_dict['weights'] = weights
                ktable_dict['ngauss'] = ngauss

            ktable_dict['kcoeff'][mol_val] = ktable['kcoeff'][:,Tmin_idx:Tmax_idx,:,:] / 10000. # from cm^-2 to m^-2

        logging.info('Loaded temperatures in ktables: %s' % ktable_dict['t'] )
        logging.info('Loaded pressures in ktables: %s ' % ktable_dict['p'])
        logging.info('The wavenumber range of the ktables is %.2f - %.2f with %i bins' % (np.min(bin_edges),
                                                                                np.max(bin_edges), len(bin_centers)))

        # set 'native' wavenumber grid
        self.int_wngrid_native = bin_centers # full wavenumber range
        with np.errstate(divide='ignore'):
            self.int_wlgrid_native = 10000./self.int_wngrid_native
        self.int_nwngrid_native = len(bin_centers)
        self.int_nwlgrid_native = len(bin_centers)
        self.ktable_dict = ktable_dict

    def load_sigma_highres_dict(self):

        logging.info('Load high resolution cross sections.')

        if self.params.atm_active_gases[0] == 'FILE':
            molecules = self.active_gases_file
        else:
            molecules = self.params.atm_active_gases

        sigma_dict = {}
        sigma_dict['xsecarr'] = {}

        temp_list_check = []
        all_molpaths = []

        logging.info('Check that the same temperatures for all molecules are available.')

        for mol_idx, mol_val in enumerate(molecules):

            # select xsec for given molecule and all temperatures available
            molpaths = insensitive_glob(os.path.join(self.params.in_xsec_path, '%s_T*' % mol_val))
            if len(molpaths) == 0:
                logging.error('There is no cross section for %s. Path: %s ' % (mol_val, self.params.in_xsec_path))
                exit()

            # get list of temperatures available
            temp_list = []
            for molpath_idx, molpath_val in enumerate(molpaths):
                # check filename
                tempfield = os.path.basename(molpath_val).split('.')[0].split('_')[1]
                if tempfield[0] != 'T':
                    logging.error('Filename of the cross section is not valid. Should be of the form H2O_T600.TauREx.pickle.'
                                  ' Check the filename for %s ' % os.path.basename(molpath_val))
                    exit()
                # get temperature of given file
                temp_list.append(float(tempfield[1:]))

            # restrict list of temperature to those needed
            temp_list_cut, Tmin_idx, Tmax_idx = self.get_temp_range_idx(np.sort(temp_list))

            # check that the restricted list of temperature is the same for all molecules
            if mol_idx > 0:
                if temp_list_cut != temp_list_check:
                    logging.error('The list of temperatures for all molecules should be the same. %s and %s don\'t share the'
                                  ' same temperatures for the range needed to compute the spectrum. ' %
                                  (molecules[mol_idx], molecules[mol_idx-1]))
                    exit()
                temp_list_check = temp_list_cut

            # loop through temperature list and load all xsec filenames (for all molecules)  into a list
            # used to get the largest file, hence the finest wavenumber grid
            # NOTE that here we assume that the largest file has the finest grid!

            for temp in temp_list_cut:
                molpath_val = insensitive_glob(os.path.join(self.params.in_xsec_path, '%s_T%i*' % (mol_val, int(temp))))[0]
                all_molpaths.append(molpath_val)

        logging.info('Beginning loading of cross sections, with interpolation to finest grid. Might take a while...')

        # get largest file in all_molpaths, and get the wavenumber grid. All other xsec will be reinterpolated to this grid
        largest_file = heapq.nlargest(1, all_molpaths, key=os.path.getsize)[0]
        largest_xsec_load = pickle.load(open(largest_file, 'rb'), encoding='latin1')
        wngrid = largest_xsec_load['wno']
        press_list = np.asarray(largest_xsec_load['p'])

        sigma_array = np.zeros((len(press_list), len(temp_list_cut), len(wngrid)))

        # loop again through all molecules and temperatures
        sigma_dict = {}
        sigma_dict['xsecarr'] = {}
        sigma_dict['t'] = np.asarray(temp_list_cut).astype(float)
        sigma_dict['p'] = press_list.astype(float)
        sigma_dict['wno'] = wngrid

        for mol_idx, mol_val in enumerate(molecules):
            logging.info('Doing %s...' % mol_val)
            for temp_idx, temp_val in enumerate(temp_list_cut):
                molpath_val = insensitive_glob(os.path.join(self.params.in_xsec_path, '%s_T%i*' % (mol_val, int(temp))))[0]

                xsec_load = pickle.load(open(molpath_val, 'rb'), encoding='latin1')

                if np.shape(xsec_load['xsecarr'])[1] != 1:
                    logging.error('This cross section has more than one temperature. For high resolution cross section'
                                  ' only one temperature per file is required. Filename: %s. ' % os.path.basename(molpath_val))

                # check that all xsec for all molecules and temperatures have the same pressures

                if not np.array_equal(np.asarray(xsec_load['p']), press_list):
                    logging.error('The cross section %s does not share the same pressure list of %s. '
                                  ' %s' % (os.path.basename(molpath_val), os.path.basename(largest_file)))
                    exit()

                # inteprolate pressure
                for press_idx, press_val in enumerate(press_list):
                    sigma_array[press_idx, temp_idx, :] = np.interp(wngrid, xsec_load['wno'], xsec_load['xsecarr'][press_idx, 0, :])

            sigma_dict['xsecarr'][mol_val] = sigma_array / 10000.  # from cm^-2 to m^-2

        # set 'native' wavenumber grid
        self.int_wngrid_native = wngrid
        self.int_nwngrid_native = len(wngrid)
        with np.errstate(divide='ignore'):
            self.int_wlgrid_native = 10000./wngrid
        self.int_nwlgrid_native = len(wngrid)
        self.sigma_dict = sigma_dict

    def load_sigma_sampled_dict(self):

        logging.info('Load sampled cross sections.')

        # list of molecules to load
        if self.params.atm_active_gases[0] == 'FILE': # loaded from external file
            molecules = self.active_gases_file
        else:
            molecules = self.params.atm_active_gases # loaded from param file

        sigma_dict = {}
        sigma_dict['xsecarr'] = {}

        for mol_idx, mol_val in enumerate(molecules):

            # check that xsec for given molecule exists
            molpath = insensitive_glob(os.path.join(self.params.in_xsec_path, '%s_*' % mol_val))+\
                      insensitive_glob(os.path.join(self.params.in_xsec_path, '%s.*' % mol_val))
            if len(molpath) > 0:
                molpath = molpath[0]
            else:
                logging.error('There is no cross section for %s. Path: %s ' % (mol_val, self.params.in_xsec_path))
                exit()

            # load cross sections
            try:
                sigma_tmp = pickle.load(open(molpath, 'rb'), encoding='latin1')
            except:
                sigma_tmp = pickle.load(open(molpath))

            # check that the wavenumber, temperature and pressure grid are the same for all cross sections
            if mol_idx > 0:
                if np.unique(sigma_tmp['wno'] - wno) != 0:
                    logging.error('The cross section of %s has a different wavenumber grid than the one of %s' %
                                  (molecules[mol_idx], molecules[mol_idx-1]))
                    exit()
                if np.unique(sigma_tmp['t'] - t) != 0:
                    logging.error('The cross section of %s has a different temperature grid than the one of %s' %
                                  (molecules[mol_idx], molecules[mol_idx-1]))
                    exit()
                if np.unique(sigma_tmp['p'] - p) != 0:
                    logging.error('The cross section of %s has a different pressure grid than the one of %s' %
                                  (molecules[mol_idx], molecules[mol_idx-1]))
                    exit()

            t = sigma_tmp['t']
            p = sigma_tmp['p']
            wno = sigma_tmp['wno']

            if mol_idx == 0:

                # restrict temperature range
                T_list, Tmin_idx, Tmax_idx = self.get_temp_range_idx(t)
                sigma_dict['t'] = T_list.astype(float)
                sigma_dict['p'] = p.astype(float)
                sigma_dict['wno'] = sigma_tmp['wno']

            sigma_dict['xsecarr'][mol_val] = sigma_tmp['xsecarr'][:,Tmin_idx:Tmax_idx] / 10000. # from cm^-2 to m^-2
        logging.info('Temperature range: %s' % sigma_dict['t'] )
        logging.info('Pressure range: %s ' % sigma_dict['p'])

        logging.info('The wavenumber range of the cross sections is %.2f - %.2f with %i points' % (np.min(wno), np.max(wno), len(wno)))

        # set 'native' wavenumber grid
        self.int_wngrid_native = wno
        self.int_nwngrid_native = len(wno)
        with np.errstate(divide='ignore'):
            self.int_wlgrid_native = 10000./wno
        self.int_nwlgrid_native = len(wno)

        self.sigma_dict = sigma_dict

    def load_sigma_rayleigh_dict(self):

        # precalculate rayleigh scattering cross sections

        logging.info('Compute Rayleigh scattering cross sections')

        sigma_rayleigh_dict = {}

        if self.params.atm_active_gases[0] == 'FILE':
            molecules = self.active_gases_file + ['H2', 'HE', 'N2']
        else:
            molecules = self.params.atm_active_gases + self.params.atm_inactive_gases

        for gasname in molecules:

            gasname = gasname.upper()

            # get the refractive index. Formulae taken from Allen Astrophysical Quantities if not otherwise specified
            n_formula = True # assume we have a formula for the refractive index of gasname
            king = 1 # King correction factor
            ns = 0   # refractive index
            wn = self.int_wngrid_native # wavenumber in cm^-1

            with np.errstate(divide='ignore'):
                wltmp = 10000./wn
            if gasname == 'HE':
                 # C. R. Mansfield and E. R. Peck. Dispersion of helium, J. Opt. Soc. Am. 59, 199-203 (1969)
                ns = 1 + 0.01470091/(423.98-(wltmp)**-2)
                # king is one for He
            elif gasname == 'H2':
                ns = 1 + 13.58e-5 * (1. + 7.52e-3 / (wltmp)**2)
                delta = 0.035 # from Morgan old code..
                king = (6.+3.*delta)/(6.-7.*delta) # Bates (1984)
            elif gasname == 'N2':
                ns = 1. + (6498.2 + (307.43305e12)/(14.4e9 - wn**2))*1.e-8 # Peck and Khanna
                king = 1.034+3.17e-12*wn**2 # Bates
            elif gasname == 'O2':
                #  J. Zhang, Z. H. Lu, and L. J. Wang. Appl. Opt. 47, 3143-3151 (2008)
                ns = 1 + 1.181494e-4 + 9.708931e-3/(75.4-(wltmp)**-2)
                king = 1.096
            elif gasname == 'CO2':
                #A. Bideau-Mehu, Y. Guern, R. Abjean and A. Johannin-Gilles. Opt. Commun. 9, 432-434 (1973)
                ns = 1 + 6.991e-2/(166.175-(wltmp)**-2)+1.44720e-3/(79.609-(wltmp)**-2)+6.42941e-5/(56.3064-(wltmp)**-2)\
                     +5.21306e-5/(46.0196-(wltmp)**-2)+1.46847e-6/(0.0584738-(wltmp)**-2)
                king = 1.1364 #  Sneep & Ubachs 2005
            elif gasname == 'CH4':
                ns = 1 + 1.e-8*(46662. + 4.02*1.e-6*(1/((wltmp)*1.e-4))**2)
            elif gasname == 'CO':
                ns = 1 + 32.7e-5 * (1. + 8.1e-3 / (wltmp)**2)
                king = 1.016  #  Sneep & Ubachs 2005
            elif gasname == 'NH3':
                ns = 1 + 37.0e-5 * (1. + 12.0e-3 / (wltmp)**2)
            elif gasname == 'H2O':
                # P. E. Ciddor. Appl. Optics 35, 1566-1573 (1996)
                ns_air = (1 + (0.05792105/(238.0185 - (wltmp)**-2) + 0.00167917/(57.362-(wltmp)**-2)))
                ns = 0.85 * (ns_air - 1.) + 1  # ns = 0.85 r(air) (Edlen 1966)
                delta = 0.17 # Marshall & Smith 1990
                king = (6.+3.*delta)/(6.-7.*delta)
            else:
                # this sets sigma_R to zero for all other gases
                n_formula = False
                logging.warning('There is no formula for the refractive index of %s. '
                                'Cannot compute the cross section' % gasname)


            if n_formula: # only if the refractive index was computed
                Ns = 2.6867805e25 # in m^-3
                with np.errstate(divide='ignore'):
                    sigma =  ( 24.*np.pi**3/ (Ns**2) )*(((ns**2-1.)/(ns**2+2.))**2) * king / ((10000./wn) * 1.e-6)**4

                # override H2 and He sigma with formulae from M Line
                with np.errstate(divide='ignore'):
                    wave = (1/wn)*1E8
                if gasname == 'H2':
                    sigma = ((8.14E-13)*(wave**(-4.))*(1+(1.572E6)*(wave**(-2.))+(1.981E12)*(wave**(-4.))))*1E-4
                if gasname == 'HE':
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

            try:
                sigma_tmp = pickle.load(open(cia_path, 'rb'), encoding='latin1') # python 3
            except:
                sigma_tmp = pickle.load(open(cia_path)) # python 2

            t = sigma_tmp['t']
            wno = sigma_tmp['wno']

            # check cia wavenumber boundaries
            if np.min(wno) > np.min(self.int_wngrid_native) or np.max(wno) < np.max(self.int_wngrid_native):
                logging.warning('Internal wavenumber grid overflow for CIA xsec for %s' % pair_val)
                logging.warning('Internal (native) wavenumber grid range: %f - %f' % (np.min(self.int_wngrid_native),
                                                                             np.max(self.int_wngrid_native)))
                logging.warning('CIA cross section wavenumber grid range: %f - %f' % (np.min(wno), np.max(wno)))
                logging.warning('Assume cia to be zero outside the xsec range')

            # restrict temperature range
            if self.params.mode == 'forward_model':
                T_list = t
                Tmin_idx = 0
            else:
                T_list, Tmin_idx, Tmax_idx = self.get_temp_range_idx(t)
            sigma_dict['t'] = T_list
            sigma_dict['wno'] = self.int_wngrid_native
            sigma_dict['xsecarr'][pair_val] = np.zeros((len(sigma_dict['t']), self.int_nwngrid_native))

            # reinterpolate cia xsec to molecule cross section grid and save
            for t_idx, t_val in enumerate(sigma_dict['t']):
                sigma_dict['xsecarr'][pair_val][t_idx,:] = np.interp(self.int_wngrid_native, wno,
                                                                     sigma_tmp['xsecarr'][Tmin_idx+t_idx])

            # load the sigma array in memory
            logging.info('Preload cia cross section for %s' % pair_val)

        del sigma_tmp, t, wno

        self.sigma_cia_dict = sigma_dict

    def load_mie_indices(self):
        '''
        Loading the refractive indices (real/imaginary) for the BH Mie model. 
        The input path is given by mie_path in parameter file. 
        File format: 
        #comment
        microns nreal nimag
        '''
        
        #loading file 
        mie_raw = np.loadtxt(self.params.in_mie_path,skiprows=1)
        
        #interpolating to native grid 
#         mie_interp = np.zeros((self.int_nwlgrid_native,2), dtype=np.float64, order='C')
#         mie_interp[:,0] = np.interp(self.int_wlgrid_native[::-1],mie_raw[:,0],mie_raw[:,1]) #real
#         mie_interp[:,1] = np.interp(self.int_wlgrid_native[::-1],mie_raw[:,0],mie_raw[:,1]) #imaginary
        
        #saving to memory
        species_name = self.params.in_mie_path.split('/')[-1].split('.')[0]
        logging.info('Preloading Mie refractive indices for %s' % species_name)
        self.mie_indices = mie_raw
        self.mie_species = species_name

        

    def load_wavenumber_grid(self):

        # restrict the native grid (i.e. the spectral grid of the loaded xsec/ktables) to various grids

        logging.info('Create wavenumber grids')

        # native. This grid is set in load_sigma_sampled_dict, load_sigma_highres_dict or load_ktables_dict
        # It is the native spectral range of the loaded cross sections or ktables. All other grids are subgrids of the native grid
        # It is used in create_spectrum if General->manual_waverange = False
        # It is also used in the output class General->manual_waverange = False to create the extended fitted spectrum.
        # The extended fitted spectrum is a spectrum extending beyond the range of obs_spectrum.
        logging.info('`native` wavenumber grid: %.3f - %.3f' % (np.min(self.int_wngrid_native), np.max(self.int_wngrid_native)))

        # obs_spectrum. Used during fitting. This is the optimal grid to fit the spectrum.
        if self.params.mode == 'retrieval':
            numin = 10000/(self.obs_wlgrid[0] + self.obs_binwidths[0]/2.)
            numax = 10000/(self.obs_wlgrid[-1] - self.obs_binwidths[-1]/2.)
            logging.info('`obs_spectrum` wavenumber grid: %.3f - %.3f'% (numin, numax))
            idx_min, idx_max = self.get_native_grid_range_idx(numin, numax, 'obs_spectrum')
            self.int_wngrid_obs_idxmin = idx_min
            self.int_wngrid_obs_idxmax = idx_max
            self.int_wngrid_obs = self.int_wngrid_native[self.int_wngrid_obs_idxmin:self.int_wngrid_obs_idxmax]
            self.int_nwngrid_obs = len(self.int_wngrid_obs)
            self.int_wlgrid_obs = 10000./self.int_wngrid_obs
            self.int_nwlgrid_obs = len(self.int_wlgrid_obs)

            self.intsp_bingrid, self.intsp_bingrididx = get_specbingrid(self.obs_wlgrid,
                                                                        self.int_wlgrid_obs,
                                                                        self.obs_binwidths)

        # manual. Used only if General->manual_waverange = True.
        # Used in create_spectrum (forward model creation), and in output class to create the extended fitted spectrum.
        # (The extended fitted spectrum is a spectrum extending beyond the range of obs_spectrum.)
        if self.params.gen_manual_waverange:
            numin = 10000/(self.params.gen_wavemax)
            numax = 10000/(self.params.gen_wavemin)
            logging.info('`manual` wavenumber grid:'
                         ' %.3f - %.3f'% (numin, numax))
            idx_min, idx_max = self.get_native_grid_range_idx(numin, numax, 'manual')
            self.int_wngrid_manual_idxmin = idx_min
            self.int_wngrid_manual_idxmax = idx_max
            self.int_wngrid_manual = self.int_wngrid_native[self.int_wngrid_manual_idxmin:self.int_wngrid_manual_idxmax]
            self.int_nwngrid_manual = len(self.int_wngrid_manual)
            self.int_wlgrid_manual = 10000./self.int_wngrid_manual
            self.int_nwlgrid_manual = len(self.int_wlgrid_manual)

        # extended (obs_spectrum extended by 10% on each side)
        # Used only after fitting, in output class. It recomputes the fitted spectrum on an extended grid with respect
        # to obs_spectrum. Also used to create 1 sigma spectrum spread.
        if self.params.mode == 'retrieval':
            numin = self.int_wngrid_obs[0] - 0.1 * self.int_wngrid_obs[0]
            numax = self.int_wngrid_obs[-1] + 0.1 * self.int_wngrid_obs[-1]
            logging.info('`extended` wavenumber grid: %.3f - %.3f' % (numin, numax))
            idx_min, idx_max = self.get_native_grid_range_idx(numin, numax, 'extended')
            self.int_wngrid_extended_idxmin = idx_min
            self.int_wngrid_extended_idxmax = idx_max
            self.int_wngrid_extended = self.int_wngrid_native[self.int_wngrid_extended_idxmin:self.int_wngrid_extended_idxmax]
            self.int_nwngrid_extended = len(self.int_wngrid_extended)
            self.int_wlgrid_extended = 10000./self.int_wngrid_extended
            self.int_nwlgrid_extended = len(self.int_wlgrid_extended)

        # # wavenumber grid limits of internal model
        # if self.params.gen_manual_waverange or not isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
        #     logging.info('Wavenumber grid limits defined in manual_waverange: %.3f - %.3f' % (10000./self.params.gen_wavemax,
        #                                                                                      10000./self.params.gen_wavemin))
        #     # limits defined by a manual wavelength range in micron in param file
        #     lambdamax = self.params.gen_wavemax
        #     lambdamin = self.params.gen_wavemin
        # else:
        #     logging.info('Wavenumber grid limits by the observed spectrum: %.3f - %.3f' % (10000./self.params.gen_wavemax,
        #                                                                                   10000./self.params.gen_wavemin))
        #     lambdamax = self.obs_wlgrid[0] + self.obs_binwidths[0]/2.
        #     lambdamin = self.obs_wlgrid[-1] - self.obs_binwidths[-1]/2.
        #
        # # convert to wavenumbers
        # numin = 10000./lambdamax
        # numax = 10000./lambdamin
        #
        # # find numin / numax closest to the cross section wavenumber grid (approx numin for defect, and numax for excess)
        # idx_min = np.argmin(np.abs(self.int_wngrid_native-numin))
        # if numin - self.int_wngrid_native[idx_min] < 0:
        #     idx_min -= 1
        #
        # idx_max = np.argmin(np.abs(self.int_wngrid_native-numax))
        # if numax - self.int_wngrid_native[idx_max] > 0:
        #     idx_max += 1
        #
        # if idx_min < 0 or idx_max < 0:
        #     logging.error('There is a problem with the internal wavenumber grid. Maybe  the spectrum is outside the '
        #                   'spectral range of the cross sections / ktables?')
        #     exit()
        #
        # self.int_wngrid_obs_idxmin = idx_min
        # self.int_wngrid_obs_idxmax = idx_max
        # self.int_wngrid_obs = self.int_wngrid_native[self.int_wngrid_obs_idxmin:self.int_wngrid_obs_idxmax]
        # self.int_nwngrid_obs = len(self.int_wngrid_obs)
        #
        # logging.info('Internal wavenumber grid is %.2f - %.2f with %i points' %
        #              (self.int_wngrid_native[self.int_wngrid_obs_idxmin],
        #               self.int_wngrid_native[self.int_wngrid_obs_idxmax-1],
        #               self.int_nwngrid_obs))
        #
        # # convert wavenumber grid to wavelenght grid
        # self.int_wlgrid_obs = 10000./self.int_wngrid_obs
        # self.int_nwlgrid_obs = len(self.int_wlgrid_obs)
        #
        # if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
        #     # calculate spectral binning grid in wavelength space
        #
        #     self.intsp_bingrid, self.intsp_bingrididx = get_specbingrid(self.obs_wlgrid,
        #                                                                 self.int_wlgrid_obs,
        #                                                                 self.obs_binwidths)
        #     self.intsp_nbingrid = len(self.obs_wlgrid)
        #     self.intsp_bingrid_full, self.intsp_bingrididx_full = get_specbingrid(self.obs_wlgrid,
        #                                                                           self.int_wlgrid_native,
        #                                                                           self.obs_binwidths)
        #     self.intsp_nbingrid_full = len(self.obs_wlgrid)


    def get_native_grid_range_idx(self, numin, numax, gridname):

        # restrict native wavenumber range to numin, numax. Return indexes of new range boundaries in native grid
        idx_min = np.argmin(np.abs(self.int_wngrid_native-numin))
        if numin - self.int_wngrid_native[idx_min] < 0:
            idx_min -= 1
        idx_max = np.argmin(np.abs(self.int_wngrid_native-numax))
        if numax - self.int_wngrid_native[idx_max] > 0:
            idx_max += 1

        if idx_min < 0 or idx_max < 0:
            logging.error('There is a problem with the %s wavenumber grid. Maybe  the range is outside the '
                          'spectral range of the `native` grid (i.e. xsec, ktables spectral range) ?'  % gridname)
            exit()

        return idx_min, idx_max

    def get_temp_range_idx(self, t):

        # restrict temperature range in input list (t)
        # return new list of temperatures, and idx of temperature boundaries of original list

        Tmin = 0
        Tmax = 5000

        if not self.params.in_custom_temp_range in ['False', 'None', None]:
            # range defined in parameter file
            Tmin =  self.params.in_custom_temp_range[0]
            Tmax =  self.params.in_custom_temp_range[1]
        else:
            # get range from other conditions (todo improve)
            # if self.params.ven_load:
            #     Tmax = np.max(self.ven_temperature)
            #     Tmin = np.min(self.ven_temperature)
            if self.params.downhill_run or self.params.mcmc_run or self.params.nest_run:
                if self.params.atm_tp_type == 'isothermal':
                    Tmax = self.params.fit_tp_iso_bounds[1]
                    Tmin = self.params.fit_tp_iso_bounds[0]
            elif self.params.atm_tp_type == 'isothermal':
                Tmin = Tmax = self.params.atm_tp_iso_temp
            if Tmax > np.max(t) or Tmin < np.min(t):
                logging.warning('The atmospheric temperature profile falls outside the temperature '
                                'range of the cross sections')
                logging.warning('Internal temperature range: %i - %i' % (Tmin, Tmax))
                logging.warning('Cross section temperature range: %i - %i' % (np.min(t), np.max(t)))

        # get idx of the temperature boundaries in the cross sections
        # always get the max T availble in the xsec that is lower than Tmin and the min T that is higher than Tmax

        if Tmax > np.max(t):
            Tmax_idx = len(t)
        else:
            Tmax_diff = Tmax - np.asarray(t) # t is sigma_tmp['t']
            Tmax_idx = len(Tmax_diff) - np.argmax(Tmax_diff[::-1][Tmax_diff<=0])

        if Tmin < np.min(t):
            Tmin_idx = 0
        else:
            Tmin_diff = Tmin - np.asarray(t) # t is sigma_tmp['t'], the list of temperature
            Tmin_idx = np.argmin(Tmin_diff[Tmin_diff>=0])

        if Tmin > np.max(t) and Tmax > np.max(t):
            Tmax_idx = len(t)
            T_list = np.asarray([t[-1]])
        else:
            T_list = t[Tmin_idx:Tmax_idx]

        return T_list, Tmin_idx, Tmax_idx


    def load_star_SED(self):

        # reading in phoenix spectra from folder specified in parameter file
        all_files = insensitive_glob(os.path.join(self.params.in_star_path, '*.fmt'))

        temperatures = []
        for filenm in all_files:
            temp = np.float(os.path.basename(filenm).split('-')[0][3:])
            temperatures.append(temp)
        temperatures = np.sort(temperatures)

        # reading in stellar file
        if self.params.star_use_blackbody or (self.params.star_temp > max(temperatures) or
                                                      self.params.star_temp < min(temperatures)):
            if self.params.verbose and not self.params.star_use_blackbody:
                logging.warning('Stellar temp. in .par file exceeds range %.1f - %.1f K. ' % (min(temperatures),
                                                                                              max(temperatures)))

            logging.info('Using black-body approximation for stellar spectrum')

            self.star_blackbody = True
            SED = black_body(self.int_wngrid_native, self.params.star_temp)
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = find_nearest(temperatures, self.params.star_temp)
            self.star_blackbody = False

            for file in all_files: #this search is explicit due to compatibility issues with Mac and Linux sorting
                if np.int(file.split('/')[-1][3:8]) == np.int(tmpselect):
                    self.SED_filename = file

            #reading in correct file and interpolating it onto self.int_wngrid_obs
            SED_raw = np.loadtxt(self.SED_filename, dtype='float', comments='#')
            SED_raw[:,1] *= 10.0 # #converting from ergs to SI
            SED = np.interp(self.int_wlgrid_native, SED_raw[:,0], SED_raw[:,1])

        self.star_sed_native = SED

    def load_active_gases_file(self):

        logging.info('Load active gases mixing ratios from external file: %s' % self.params.atm_active_gases_file)

        # read list of molecules (first line)
        with open(self.params.atm_active_gases_file, 'r') as f:
            first_line = f.readline()
            molecules = ' '.join(first_line.split()).split(' ')

        # read mixing ratio profiles (from second row)
        if os.path.isfile(self.params.atm_active_gases_file):
            mixratio_profiles_file = np.loadtxt(self.params.atm_active_gases_file, skiprows=1)

        # first column is pressure
        self.mixratio_pressure_file = mixratio_profiles_file[:,0]
        nrows = len(self.mixratio_pressure_file) # number of pressure levels

        # set inactive gases. Allow H2, HE and N2. If one of these is not present, just set to zero
        self.inactive_gases_file = ['H2', 'HE', 'N2']
        self.inactive_mixratio_profile_file = np.zeros((nrows, 3))
        count_inactive = 0
        for inactive_mol_idx, inactive_mol_val in enumerate(['H2', 'HE', 'N2']):
            for mol_idx, mol_val in enumerate(molecules):
                if mol.upper() == inactive_mol:
                    count_inactive += 1
                    self.inactive_mixratio_profile_file[:,inactive_mol_idx] = mixratio_profiles_file[:,mol_idx+1] # plus 1 as first column is pressure

        # Set active gases (get everything, excluding H2, HE, N2
        self.active_gases_file = []
        self.active_mixratio_profile_file = np.zeros((nrows, len(molecules) - count_inactive))
        idx = 0
        for mol_idx, mol_val in enumerate(molecules):
            if not mol_val.upper() in ['H2', 'HE', 'N2']:

                # check that we have cross section or ktable for this molecule
                if self.opacity_method[:4] == 'xsec':
                    opacity_path = self.params.in_xsec_path
                elif self.opacity_method == 'ktables':
                    opacity_path = self.params.in_ktab_path

                # check that xsec for given molecule exists
                molpath = insensitive_glob(os.path.join(opacity_path, '%s_*' % mol_val))+\
                          insensitive_glob(os.path.join(opacity_path, '%s.*' % mol_val))
                if len(molpath) <= 0: # we don't have xsec/ktable for this moluecule
                    excluded_molecules.append(mol_val)

                else:
                    self.active_gases_file.append(mol_val)
                    self.active_mixratio_profile_file[:,idx] = mixratio_profiles_file[:,mol_idx+1]
                    idx += 1

        if len(excluded_molecules) > 1:
            logging.warning('There are some excluded molecules from the input mixing ratios:', excluded_molecules)

    def load_tp_profile_file(self):

        logging.info('Load temperature profile from external file: %s' % self.params.atm_tp_file)

        # note first column is pressure, secondo column is temperature
        if os.path.isfile(self.params.atm_tp_file):
            self.tp_profile_file = np.loadtxt(self.params.atm_tp_file)