'''
   TauREx v2 - Development version

    Data class
    
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
            self.obs_spectrum = self.readfile(params.in_spectrum_file)
        else:
            # no input spectrum provided. If running only forward model, input spectrum can be omitted
            logging.info('No input spectrum provided')
            self.obs_spectrum = False
        # observed input spectrum
        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # set observed spectrum specific variables (only if spectrum is provided)
            self.obs_wlgrid = self.obs_spectrum[:,0] # wavegrid in micron
            self.obs_wngrid = 10000./self.obs_spectrum[:,0]
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
        #approximate numin and numax to closest gridding point
        # todo NOT WORKING PROPERLY for in_xsec_dnu < 1 !!!
        numin = round_base(numin, self.params.in_xsec_dnu)
        numax = round_base(numax, self.params.in_xsec_dnu)
        # create wavenumber grid of internal model using numin, numax and delta wavenumber provided in param file
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

        #reading in atmospheric profile file
        if self.params.in_use_ATMfile:
            logging.info('Using .atm file to generate atmospheric TP profile')
            self.pta, self.X = self.readATMfile() # pta = pressure, temp, alt; X = mixing ratios of molecules
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])
        else:
            self.nlayers = self.params.tp_atm_levels
            self.ngas = len(self.params.planet_molec)

        # Rayleigh scattering
        if MPIrank == 0: #
            # compute rayleigh scattering cross sections and save them to a file
            self.build_rayleigh_sigma()
        self.sigma_R = self.get_rayleigh_sigma() # load Rayleigh scattering cross sections

        # read collision induced absorption file
        if params.in_include_cia:
            self.cia = self.readfile(self.params.in_cia_file, interpolate=True)
        
        #reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.F_star = self.get_star_SED() #@todo there is most certainly a bug there. think units are ergs/s/cm^2 at the moment

        logging.info('Loading cross sections into arrays')

        molecules = self.params.planet_molec
        sigmadb_all = [pickle.load(open(os.path.join(self.params.in_xsec_path, '%s.db' % molecule))) for molecule in molecules]
        self.sigma_nmol = len(molecules)
        self.sigma_wno = np.asarray([sigmadb['wno'] for sigmadb in sigmadb_all])
        self.sigma_np = np.asarray([len(sigmadb['p']) for sigmadb in sigmadb_all])
        self.sigma_nt = np.asarray([len(sigmadb['t']) for sigmadb in sigmadb_all])
        self.sigma_p = np.asarray([sigmadb['p'] for sigmadb in sigmadb_all])
        self.sigma_t = np.asarray([sigmadb['t'] for sigmadb in sigmadb_all])
        self.sigma_array = []
        int_xsec_sigma = np.zeros(self.int_nwngrid) # assume sigma = 0 for regions of the internal wn grid not covered by the xsec
        for idx_m, val_molecule in enumerate(molecules):
            logging.info('Reading %s: %i pressures, %i temperatures between %.2f and %.2f cm-1' % (val_molecule, self.sigma_np[idx_m], self.sigma_nt[idx_m],
                                                                                               np.min(self.sigma_wno[idx_m]), np.max(self.sigma_wno[idx_m])))
            ext_xsec_wno = self.sigma_wno[idx_m]
            for idx_p in xrange(self.sigma_np[idx_m]):
                for idx_t in xrange(self.sigma_nt[idx_m]):
                    # adapt input cross section wavenumber grid to internal grid
                    ext_xsec_sigma = np.asarray(sigmadb_all[idx_m]['xsecarr'][idx_p, idx_t, :])
                    idmin_xsec = 0
                    idmin_sigdic = 0
                    if np.min(ext_xsec_wno) < np.min(self.int_wngrid):
                        idmin_xsec = np.where(np.abs(ext_xsec_wno - np.min(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]
                    else:
                        idmin_sigdic = np.where(np.abs(self.int_wngrid - np.min(ext_xsec_wno)) < self.params.in_xsec_dnu)[0][0]
                    idmax_xsec = len(ext_xsec_wno) - 1
                    idmax_sigdic = self.int_nwngrid - 1
                    if np.max(ext_xsec_wno) > np.max(self.int_wngrid):
                        idmax_xsec = np.where(np.abs(ext_xsec_wno - np.max(self.int_wngrid)) < self.params.in_xsec_dnu)[0][0]
                    else:
                        idmax_sigdic = np.where(np.abs(self.int_wngrid - np.max(ext_xsec_wno)) < self.params.in_xsec_dnu)[0][0]
                    int_xsec_sigma =  (ext_xsec_sigma[idmin_xsec:idmax_xsec] * 1e-4).tolist() # convert cm^2 to m^2
                    self.sigma_array = self.sigma_array + int_xsec_sigma
        self.sigma_array = np.asarray(self.sigma_array)
            
        logging.info('Data object initialised')

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

    def readfile(self, name, interpolate=False):
        try:
            out = loadtxt(name)
        except ValueError:
            out = loadtxt(name, delimiter=',')
        if len(out[:,0]) < len(out[0,:]):
            out = transpose(out)
        #sorting data along ascending first column
        out = out[argsort(out[:,0]),:]
        if interpolate:
            # inteprolate to wavenumber grid
            interpflux = interp(self.int_wngrid, out[:,0], out[:,1])
            out = transpose(vstack((self.int_wngrid, interpflux)))
        return out

    def set_ABSfile(self,path=None,filelist=None,interpolate = False,temperature=None):
        # @todo  Preselector stuff?!?!? CHECK!
        #manually overwrites absorption coefficients from new file
        #input path needs to be given and list of filename strings
        if path is None:
            extpath = self.params.in_abs_path
        if filelist is None:
            raise IOError('No input ABS file specified')
        if temperature is None:
            temperature = self.params.planet_temp/Users/marco/Dropbox/repos/TauREx/library/pathintegral.cpp
        sigma_array,self.int_wngrid = self.readABSfiles(extpath=path,
                                                      extfilelist=filelist,
                                                      interpolate2grid=interpolate,
                                                      outputwavegrid=True)
        self.sigma_dict = {}
        self.sigma_dict['tempgrid'] = [int(temperature)]
        self.sigma_dict[int(temperature)] = sigma_array
        self.int_nwngrid = len(self.int_wngrid)

    def readATMfile(self):
        
        try:
            out = np.loadtxt(self.params.in_atm_file)
        except ValueError:
            out = np.loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
        out[:,2] *= 1000. #converting from km to m
        out = out[np.argsort(out[:,2]),:]
        return out[:,0:3],np.transpose(out[:,3:])

    def build_rayleigh_sigma(self):

        # write rayleigh scattering cross sections to files

        logging.info('Write rayleigh scattering cross sections to files')

        gases = self.params.planet_molec + self.params.planet_inactive_gases

        for gasname in gases:

            filename = os.path.join('Input', 'rayleigh', '%s.rayleigh' % gasname)

            if not os.path.isfile(filename):

                # get the refractive index of a given gas gasname for 1000-30000 wn grid
                # Formulae taken from Allen Astrophysical Quantities if not otherwise specified

                n_formula = True # assume we have a formula for the refractive index of gasname
                king = 1 # King correction factor
                ns = 0   # refractive index
                wn = np.arange(1, 30000, 1) # wavenumber in cm^-1

                if gasname == 'He':
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
                    sigma_R =  ( 24.*np.pi**3/ (Ns**2) ) * (((ns**2 - 1.)/(ns**2 + 2.))**2) * king / ((10000./wn) * 1.e-6)**4
                    logging.info('Saving cross section for %s to %s' % (gasname, filename))
                    np.savetxt(filename, np.vstack((wn, sigma_R)).transpose())

    def get_rayleigh_sigma(self):

        logging.info('Loading rayleigh scattering cross section from files')
        gases = self.params.planet_molec + self.params.planet_inactive_gases
        sigma_R = {}
        for gasname in gases:
            filename = os.path.join('Input', 'rayleigh', '%s.rayleigh' % gasname)
            if os.path.isfile(filename):
                data = np.loadtxt(os.path.join('Input', 'rayleigh', '%s.rayleigh' % gasname))
                sigma = interp1d(data[:,0], data[:,1])
                logging.info('Reading Rayleigh scattering cross section for gas %s' % gasname)
                sigma_R[gasname] = sigma
            else:
                logging.warning('Cannot find Rayleigh scattering cross section for gas %s. Set to 0' % gasname)
        return sigma_R

    def get_molecular_weight(self, gasname):

        if gasname == 'He':
            mu = 4.
        elif gasname == 'H2':
            mu = 2.
        elif gasname == 'N2':
            mu = 28.
        elif gasname == 'O2':
            mu = 32.
        elif gasname == 'CO2':
            mu = 44.
        elif gasname == 'CH4':
            mu = 16.
        elif gasname == 'CO':
            mu = 28.
        elif gasname == 'NH3':
            mu = 17.
        elif gasname == 'H2O':
            mu = 18.
        elif gasname == 'C2H2':
            mu = 26.04
        elif gasname == 'H2S':
            mu = 34.0809
        else:
            mu = 0

        return mu * AMU