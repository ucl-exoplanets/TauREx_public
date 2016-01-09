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
from scipy import interpolate

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

    def __init__(self, params, spectrum=None):

        logging.info('Initialising data object')

        #checking taurex license
        license_manager().run()
        
        self.params = params

        #reading in spectrum data to be fitted
        if isinstance(spectrum, (np.ndarray, np.generic)):
            # read spectrum from argument
            logging.info('Reading spectrum from argument')
            self.obs_spectrum = spectrum
        elif params.in_spectrum_file:
            # read spectrum from file
            logging.info('Reading spectrum from file')
            self.obs_spectrum = np.loadtxt(params.in_spectrum_file)
        else:
            # no input spectrum provided. If running only forward model, input spectrum can be omitted
            logging.info('No input spectrum provided')
            self.obs_spectrum = False

        # observed input spectrum
        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # set observed spectrum specific variables (only if spectrum is provided)
            self.obs_wlgrid = self.obs_spectrum[:,0] # grid in micron
            self.obs_wngrid = 10000./self.obs_spectrum[:,0] # grid in wavenumbers
            self.obs_nwlgrid = len(self.obs_spectrum[:,0]) # number of datapoints in spectrum
            self.obs_binwidths = self.obs_spectrum[:,3]   if shape(self.obs_spectrum)[1] == 4 else None # bin widths

        logging.info('Create wavenumber grid of internal model')
        
        # wavenumber grid of internal model
        if params.gen_manual_waverange or not isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # limits defined by a manual wavelength range in micron in param file
            lambdamax = self.params.gen_wavemax
            lambdamin = self.params.gen_wavemin
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

        #approximate numin and numax to closest xsec gridding point
        numin = round_base(numin, self.params.in_xsec_dnu)
        numax = round_base(numax, self.params.in_xsec_dnu)

        # create wavenumber grid of internal model using numin, numax and delta wavenumber provided in param file
        # NB: the xsec grid needs to be sampled in multiples of self.params.in_xsec_dnu
        self.int_wngrid = np.arange(numin, numax, self.params.in_xsec_dnu)
        self.int_dwngrid = np.diff(self.int_wngrid)
        self.int_nwngrid = len(self.int_wngrid)
        logging.info('Internal wavenumber grid is %.2f - %.2f in steps of %.2f, resulting in %i points' % (numin, numax, self.params.in_xsec_dnu, self.int_nwngrid))

        # convert wavenumber grid to wavelenght grid
        self.int_wlgrid = 10000./self.int_wngrid
        self.int_dwlgrid = np.diff(self.int_wlgrid)
        self.int_nwlgrid = len(self.int_wlgrid)
        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # calculate spectral binning grid in wavelength space
            self.intsp_bingrid, self.intsp_bingrididx = get_specbingrid(self.obs_wlgrid, self.int_wlgrid, self.obs_binwidths)
            self.intsp_nbingrid = len(self.obs_wlgrid)

        # #reading in atmospheric profile file

        # reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.F_star = self.get_star_SED() # todo there is most certainly a bug there. think units are ergs/s/cm^2 at the moment

        # Preload the cross sections for the given wavenumber range and temperature range
        self.sigma_dict = self.load_sigma_dict()

        # Preload rayleigh cross secitons
        self.sigma_rayleigh_dict = self.load_sigma_rayleigh_dict()

        # Preload CIA cross sections
        self.sigma_cia_dict = self.load_sigma_cia_dict()

        logging.info('Data object initialised')


    def load_sigma_dict(self):

        # preload the absorption cross sections for the input molecules

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
                    if np.min(wno) > np.min(self.int_wngrid) or np.max(wno) < np.max(self.int_wngrid):
                        logging.error('Internal wavenumber grid overflow')
                        logging.error('Internal wavenumber grid range: %f - %f' % (np.min(self.int_wngrid), np.max(self.int_wngrid)))
                        logging.error('Cross section wavenumber grid range: %f - %f' % (np.min(wno), np.max(wno)))
                        exit()
                    if np.unique(np.diff(wno)) != self.params.in_xsec_dnu:
                        logging.error('The resolution of the internal model (%f wno) is different from '
                                      'the resolution of the cross sections (%f wno)' % (self.params.in_xsec_dnu, np.unique(np.diff(wno))))
                        exit()

                    # restrict wavenumber range
                    wno_min_idx = np.where(np.abs(wno - np.min(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]
                    wno_max_idx = np.where(np.abs(wno - np.max(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]+1

                    # restrict temperature range
                    Tmax = Tmin = None
                    if self.params.downhill_run or self.params.mcmc_run or self.params.nest_run:
                        Tmax = self.params.fit_T_bounds[1]
                        Tmin = self.params.fit_T_bounds[0]
                    elif self.params.in_use_ATMfile:
                        Tmax = np.max(self.pta[:,1])
                        Tmin = np.min(self.pta[:,1])
                    else:
                        Tmin = Tmax = self.params.planet_temp

                    if Tmax > np.max(t) or Tmin < np.min(t):
                        logging.error('The atmospheric temperature profile falls outside the temperature range of the cross sections')
                        logging.error('Internal temperature range: %i - %i' % (Tmin, Tmax))
                        logging.error('Cross section temperature range: %i - %i' % (np.min(t), np.max(t)))
                        exit()
                    if Tmax and Tmin:
                        # get idx of the temperature boundaries in the cross sections
                        # always get the max T availble in the xsec that is lower than Tmin and the min T that is higher than Tmax
                        Tmin_diff = Tmin - t # t is sigma_tmp['t'], the list of temperature
                        Tmin_idx = np.argmin(Tmin_diff[Tmin_diff>=0])
                        Tmax_diff = Tmax - t # t is sigma_tmp['t']
                        Tmax_idx = len(Tmax_diff) - np.argmax(Tmax_diff[::-1][Tmax_diff<=0])
                    else:
                        Tmin_idx = 0
                        Tmax_idx = len(t) - 1

                    sigma_dict['t'] = t[Tmin_idx:Tmax_idx]
                    sigma_dict['p'] = p
                    sigma_dict['wno'] = sigma_tmp['wno'][wno_min_idx:wno_max_idx] # check: this should be identical to internal model!

                    if np.unique(sigma_dict['wno'] - self.int_wngrid) != 0:
                        logging.error('Wavenumber grid of xsec is different from wn grid of internal model')
                        exit()

                # load the sigma array in memory
                logging.info('Preload cross section for %s' % mol_val)
                sigma_dict['xsecarr'][mol_val] = sigma_tmp['xsecarr'][:,Tmin_idx:Tmax_idx,wno_min_idx:wno_max_idx] / 10000. # also convert from cm^-2 to m^-2

        del sigma_tmp, t, p, wno

        return sigma_dict

    def load_sigma_rayleigh_dict(self):

        # precalculate rayleigh scattering cross sections

        logging.info('Compute Rayleigh scattering cross sections')

        sigma_rayleigh_dict = {}

        for gasname in self.params.atm_active_gases + self.params.atm_inactive_gases:

            gasname = gasname.upper()

            # get the refractive index. Formulae taken from Allen Astrophysical Quantities if not otherwise specified
            n_formula = True # assume we have a formula for the refractive index of gasname
            king = 1 # King correction factor
            ns = 0   # refractive index
            wn = self.int_wngrid # wavenumber in cm^-1

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
                logging.info('Rayleigh scattering cross section of %s correctly computed' % (gasname))
                sigma_rayleigh_dict[gasname] = sigma

        return sigma_rayleigh_dict

    def load_sigma_cia_dict(self):

        # preload collision induced absorption cross sections

        molecules = self.params.atm_active_gases
        sigma_dict = {}
        sigma_dict['xsecarr'] = {}

        for pair_val in self.params.atm_cia_pairs:

            cia_path = os.path.join(self.params.in_cia_path, '%s.db' % pair_val.upper())

            sigma_tmp = pickle.load(open(cia_path)) # load pickled cross section array

            t = sigma_tmp['t']
            wno = sigma_tmp['wno']

            # check cia wavenumber boundaries
            if np.min(wno) > np.min(self.int_wngrid) or np.max(wno) < np.max(self.int_wngrid):
                logging.warning('Internal wavenumber grid overflow for CIA xsec for %s' % pair_val)
                logging.warning('Internal wavenumber grid range: %f - %f' % (np.min(self.int_wngrid), np.max(self.int_wngrid)))
                logging.warning('CIA cross section wavenumber grid range: %f - %f' % (np.min(wno), np.max(wno)))
                logging.warning('Assume xsec to be zero outside the cia range')

            # restrict temperature range
            Tmax = Tmin = None
            if self.params.downhill_run or self.params.mcmc_run or self.params.nest_run:
                Tmax = self.params.fit_T_bounds[1]
                Tmin = self.params.fit_T_bounds[0]
            elif self.params.in_use_ATMfile:
                Tmax = np.max(self.pta[:,1])
                Tmin = np.min(self.pta[:,1])
            else:
                Tmin = Tmax = self.params.planet_temp

            if Tmax > np.max(t) or Tmin < np.min(t):
                logging.error('The atmospheric temperature profile falls outside the temperature range of the cross sections')
                logging.error('Internal temperature range: %i - %i' % (Tmin, Tmax))
                logging.error('Cross section temperature range: %i - %i' % (np.min(t), np.max(t)))
                exit()
            if Tmax and Tmin:
                # get idx of the temperature boundaries in the cross sections
                # always get the max T availble in the xsec that is lower than Tmin and the min T that is higher than Tmax
                Tmin_diff = Tmin - t # t is sigma_tmp['t'], the list of temperature
                Tmin_idx = np.argmin(Tmin_diff[Tmin_diff>=0])
                Tmax_diff = Tmax - t # t is sigma_tmp['t']
                Tmax_idx = len(Tmax_diff) - np.argmax(Tmax_diff[::-1][Tmax_diff<=0])
            else:
                Tmin_idx = 0
                Tmax_idx = len(t) - 1

            sigma_dict['t'] = t[Tmin_idx:Tmax_idx]
            sigma_dict['wno'] = self.int_wngrid
            sigma_dict['xsecarr'][pair_val] = np.zeros((len(sigma_dict['t']), self.int_nwngrid))

            # reinterpolate cia xsec to internal grid and save
            for t_idx, t_val in enumerate(sigma_dict['t']):
                sigma_dict['xsecarr'][pair_val][t_idx,:] = np.interp(self.int_wngrid, wno, sigma_tmp['xsecarr'][Tmin_idx+t_idx])

            # load the sigma array in memory
            logging.info('Preload cia cross section for %s' % pair_val)

        del sigma_tmp, t, wno

        return sigma_dict

    def get_star_SED(self):

        #reading in phoenix spectra from folder specified in parameter file
        
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
            SED = black_body(self.int_wngrid,self.params.star_temp) #@todo bug here? not multiplied by size of star 4piRs^2
#             SED *= self.params.star_radius**2 * np.pi * 4.
#             SED *= self.params.star_radius**2 * np.pi * 4.
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = find_nearest(tmpind, self.params.star_temp)
            self.star_blackbody = False
            
            for file in fileindex: #this search is explicit due to compatibility issues with Mac and Linux sorting
                if np.int(file.split('/')[-1][3:8]) == np.int(tmpselect):
                    self.SED_filename = file

            #reading in correct file and interpolating it onto self.int_wngrid
            SED_raw = np.loadtxt(self.SED_filename, dtype='float', comments='#')
            SED_raw[:,1] *= 10.0  #converting from ergs to SI @todo move converting somewhere more sane 
#             SED_raw[:,1] *= self.params.star_radius**2 * np.pi * 4.
#             digitized = np.digitize(SED_raw[:,0],self.int_wngrid)
#             SED = np.asarray([SED_raw[digitized==i,1].mean() for i in range(0,len(self.int_wngrid))])
            SED = np.interp(self.int_wngrid, SED_raw[:,0], SED_raw[:,1])
        
#         print self.params.star_temp
#         SED = black_body(self.int_wngrid,self.params.star_temp)
        return SED


    def readATMfile(self):
        
        try:
            out = np.loadtxt(self.params.in_atm_file)
        except ValueError:
            out = np.loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
        out[:,2] *= 1000. #converting from km to m
        out = out[np.argsort(out[:,2]),:]
        return out[:,0:3],np.transpose(out[:,3:])

