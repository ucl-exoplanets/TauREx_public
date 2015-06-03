################################################
#class data
#
# Reads in all relevant data and performs pre-processing 
# such as sorting, grid-interpolations etc. 
#
# Input: -parameter object
#
#
# Output: -data object containing relevant data in dictionary 
#         format  
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013 
#
################################################


#loading libraries     
import numpy, pylab, os, glob
from base import base
from numpy import *
from pylab import *
from StringIO import StringIO
from scipy.interpolate import interp1d
from scipy import interpolate

import library_general as libgen
import library_emission as libem
import logging

#loading taurex license manager. Only loaded in data class
import license
from license import *

# some constants
KBOLTZ = 1.380648813e-23
G = 6.67384e-11
AMU   = 1.660538921e-27 #atomic mass to kg

class data(base):

    def __init__(self, params, spectrum=None):

        logging.info('Initialising data object')

        #checking taurex license
        license_manager().run()

        self.params = params

        #converting absorption cross-sectinos from cm^-1 to microns
        if params.in_convert2microns:
            libgen.convert2microns(params.in_abs_path)
            libgen.convert2microns(params.in_abs_path_P)

        #reading in spectrum data to be fitted
        if isinstance(spectrum, (np.ndarray, np.generic)):
            # read spectrum from argument
            self.spectrum = spectrum
        elif params.in_spectrum_file:
            # read spectrum from file
            self.spectrum = self.readfile(params.in_spectrum_file)
        else:
            # if running only forward model, input spectrum can be omitted
            self.spectrum = False

        if isinstance(self.spectrum, (np.ndarray, np.generic)):
            # set observed spectrum specific variables (only if spectrum is provided)
            self.nwave = len(self.spectrum[:,0])
            self.wavegrid = self.spectrum[:,0]

        #calculating wavelength grids
        if params.gen_manual_waverange:
            # manual wavelength range provided in parameter file
            if not self.params.gen_abs_wavegrid:
                # generate wavelength grid with uniform binning in log(lambda)
                self.specgrid, self.dlamb_grid = self.get_specgrid(R=self.params.gen_spec_res,
                                                                   lambda_min=self.params.gen_wavemin,
                                                                   lambda_max=self.params.gen_wavemax)
            else:
                # use the wavelength grid of the cross sections
                self.specgrid, self.dlamb_grid = self.get_specgrid_from_crosssections(lambda_min=self.params.gen_wavemin,
                                                                                      lambda_max=self.params.gen_wavemax)

        else:
            # user wavelength range from observed spectrum, but increase model wavelength range by one data spectral bin
            # in blue and red. This ensures correct model binning at the edges
            bin_up = self.wavegrid[-1]-self.wavegrid[-2]
            bin_low = self.wavegrid[1]-self.wavegrid[0]
            if not self.params.gen_abs_wavegrid:
            # generate wavelength grid with uniform binning in log(lambda)
                self.specgrid, self.dlamb_grid = self.get_specgrid(R=self.params.gen_spec_res,
                                                                   lambda_min=self.wavegrid[0]-bin_low,
                                                                   lambda_max=self.wavegrid[-1]+bin_up)
            else:
                # use the wavelength grid of the cross sections
                self.specgrid, self.dlamb_grid = self.get_specgrid_from_crosssections(lambda_min=self.wavegrid[0]-bin_low,
                                                                                      lambda_max=self.wavegrid[-1]+bin_up)

        self.nspecgrid = len(self.specgrid)

        if isinstance(self.spectrum, (np.ndarray, np.generic)):
            #calculating spectral binning grid, only if observed spectrum is provided
            self.spec_bin_grid, self.spec_bin_grid_idx = self.get_specbingrid(self.wavegrid, self.specgrid)
            self.n_spec_bin_grid= len(self.spec_bin_grid)

        #reading in atmospheric profile file
        if self.params.in_use_ATMfile:
            self.pta, self.X = self.readATMfile() # pta = pressure, temp, alt; X = mixing ratios of molecules
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])

        # Rayleigh scattering
        if self.params.in_create_sigma_rayleigh:
            # compute cross sections and save them to a file
            self.build_rayleigh_sigma()
        self.sigma_R = self.get_rayleigh_sigma() # load Rayleigh scattering cross sections from file

        # Absorption cross sections
        if self.params.gen_run_gui:
            # if running in GUI mode, the dictionary is loaded from a file
            # todo add a param to rebuild sigma_array even in GUI mode
            import pickle
            filename = 'Input/sigma_dict.pickle'
            if os.path.isfile(filename):
                # sigma dict has been presaved with the input parameters (lambda = 0.4-20, resolution = 500)
                logging.info('Loading sigma dictionary from file')
                self.sigma_dict = pickle.load(open(filename, 'rb'))
            else:
                logging.info('Creating sigma file')
                self.sigma_dict  = self.build_sigma_dic(tempstep=params.in_tempres)
                pickle.dump(self.sigma_dict, open(filename, 'wb'))

        else: # not running in GUI mode
            #reading in absorption coefficient data

            if self.params.in_use_P_broadening:
                self.sigma_dict_pres, self.sigma_templist, self.sigma_preslist  = self.build_sigma_dic_pressure()
            else:
                self.sigma_dict  = self.build_sigma_dic(tempstep=params.in_tempres)

        # #reading in other files if specified in parameter file
        if params.in_include_rad:
            self.rad = self.readfile(self.params.in_rad_file, interpolate=True)
        if params.in_include_cia:
            self.cia = self.readfile(self.params.in_cia_file, interpolate=True)
        # if params.in_include_cld:
        #     self.cld = self.readfile(self.params.in_cld_file, interpolate=True)

        #reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.F_star = self.get_star_SED() #@todo there is most certainly a bug there. think units are ergs/s/cm^2 at the moment


    #class functions

    #@profile
    def get_specgrid(self, R=5000, lambda_min=0.1, lambda_max=20.0):
        #generating wavelength grid with uniform binning in log(lambda)
        #lambda min and max are in microns, R is the spectral resolution
        #R = lambda/delta_lamdba
        specgrid = []
        delta_lambda =[]
        specgrid.append(lambda_min)
        run = True
        i=1
        while run:
            dlam= specgrid[i-1]/R
            specgrid.append(specgrid[i-1]+dlam)
            delta_lambda.append(dlam)

            if specgrid[i] >= lambda_max:
                run=False
            i+=1


        return np.asarray(specgrid),np.asarray(delta_lambda)
    
    #@profile
    def get_specbingrid(self,wavegrid, specgrid):
        #function calculating the bin boundaries for the data 
        #this is used to bin the internal spectrum to the data in fitting module
        
        bingrid =[]
        bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
        for i in range(len(wavegrid)-1):
            bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
        bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
        
        bingrid_idx = numpy.digitize(specgrid,bingrid) #getting the specgrid indexes for bins

        return bingrid, bingrid_idx 

    def get_specgrid_from_crosssections(self, lambda_min=0.1, lambda_max=20.0):
        # return the wavelength grid of the absorbion cross sections
        # assume that the cross sections have a uniform wl grid

        if self.params.in_use_P_broadening:
            path = self.params.in_abs_path_P
        else:
            path = self.params.in_abs_path

        files = glob.glob(os.path.join(path, '*.abs'))
        out = self.readfile(files[0])
        wlgrid = out[:,0]

        # get only values between lambda_min and lambda_max
        specgrid = wlgrid[np.logical_and(wlgrid>lambda_min, wlgrid<lambda_max)]
        # todo: add one bin on either side of specgrid

        dlamb_grid = np.diff(specgrid)

        return specgrid, dlamb_grid


    #@profile
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
            SED = libem.black_body(self.specgrid,self.params.star_temp) #@todo bug here? not multiplied by size of star 4piRs^2
#             SED *= self.params.star_radius**2 * np.pi * 4.

#             SED *= self.params.star_radius**2 * np.pi * 4.
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = libgen.find_nearest(tmpind, self.params.star_temp)
            self.star_blackbody = False
            
            for file in fileindex: #this search is explicit due to compatibility issues with Mac and Linux sorting
                if np.int(file.split('/')[-1][3:8]) == np.int(tmpselect):
                    self.SED_filename = file
                    print file

            #reading in correct file and interpolating it onto self.specgrid
            SED_raw = np.loadtxt(self.SED_filename, dtype='float', comments='#')
            SED_raw[:,1] *= 10.0  #converting from ergs to SI @todo move converting somewhere more sane 
#             SED_raw[:,1] *= self.params.star_radius**2 * np.pi * 4. 
            
#             digitized = np.digitize(SED_raw[:,0],self.specgrid)
#             SED = np.asarray([SED_raw[digitized==i,1].mean() for i in range(0,len(self.specgrid))])
            SED = np.interp(self.specgrid, SED_raw[:,0], SED_raw[:,1])
        
#         print self.params.star_temp
#         SED = libem.black_body(self.specgrid,self.params.star_temp)
        return SED

    #@profile
    def readfile(self, name, interpolate=False):

        #reads in data file with columns wavelength and data

        try:
            out = loadtxt(name)
        except ValueError:
            out = loadtxt(name, delimiter=',')
        if len(out[:,0]) < len(out[0,:]):
            out = transpose(out)

        #sorting data along ascending first column
        out = out[argsort(out[:,0]),:]

        #bin to specgrid
        if interpolate:
            interpflux = interp(self.specgrid, out[:,0], out[:,1])
            out = transpose(vstack((self.specgrid, interpflux)))

            # bin_grid, bin_grid_idx = self.get_specbingrid(self.specgrid, out[:,0])
            # binnedvalues = np.asarray([out[:,1][bin_grid_idx == i].mean() for i in xrange(1,len(bin_grid))])
            # binnedvalues[np.isnan(binnedvalues)] = 0
            # out = transpose(vstack((self.specgrid, binnedvalues)))

        return out

    #@profile
    def build_sigma_dic(self, tempstep=50):

        #building temperature dependent sigma_array

        mollist = self.params.planet_molec

        tempmin = []
        tempmax = []

        moldict = {}

        for molecule in mollist:

            logging.info('Load sigma array for molecule %s' % molecule)

            moldict[molecule] ={}

            absfiles, templist = libgen.find_absfiles(self.params.in_abs_path, molecule)
            tempmax.append(np.max(templist))
            tempmin.append(np.min(templist))

            sigtmp = self.readABSfiles(extfilelist=absfiles,
                                       extpath=self.params.in_abs_path,
                                       interpolate2grid=True,
                                       outputwavegrid = False)
            sigshape = np.shape(sigtmp)

            moldict[molecule]['sigma'] = sigtmp
            moldict[molecule]['templist'] = templist
            moldict[molecule]['tempmin'] = tempmin
            moldict[molecule]['tempmax'] = tempmax

        #setting up new temperature grid
        temp_totmin = np.max(tempmin)
        temp_totmax = np.min(tempmax)
        tempgrid = [i for i in range(int(temp_totmin),int(temp_totmax),int(tempstep))]
        tempgrid.append(int(temp_totmax))
        moldict['tempgrid'] = tempgrid

        #interpolating absorption cross sections to new temperature grid
        for molecule in mollist:
            interpsigma = np.zeros((len(tempgrid),sigshape[1]))

            for i in range(sigshape[1]):
                interpsigma[:,i] = np.interp(tempgrid,moldict[molecule]['templist'],moldict[molecule]['sigma'][:,i]) #linear interpolation. to be changed with hill et al method
            moldict[molecule]['interpsigma'] = interpsigma

        #building final sigma_array dictionary
        sigma_dict = {}
        sigma_dict['tempgrid'] = tempgrid
        for i in range(len(tempgrid)):
            sigma_dict[tempgrid[i]] = {}

            sigma_array = np.zeros((len(mollist),sigshape[1]))
            j = 0
            for molecule in mollist:
                sigma_array[j,:] = moldict[molecule]['interpsigma'][i,:]
                j += 1
            sigma_dict[tempgrid[i]] = sigma_array

        return sigma_dict


    def build_sigma_dic_pressure(self):

        # building temperature and pressure dependent sigma_array
        # pressstep in bar
        # assume uniform temperature pressure grids for each molecule
        # (each temperature has the same pressures, and viceversa)

        import pickle

        mollist = self.params.planet_molec
        moldict = {}

        # initialise pressure and temperature lists
        templist = None
        preslist = None

        for molecule in mollist:

            logging.info('Load sigma array for molecule %s' % molecule)
            molpath = os.path.join(self.params.in_abs_path_P, molecule)
            absfilelist, templist_tmp, preslist_tmp = libgen.find_absfiles_pressure(molpath, molecule)

            if (templist == None or templist == templist_tmp) and (preslist == None or preslist == preslist_tmp):
                templist = templist_tmp
                preslist = preslist_tmp
            else:
                logging.error('Cannot build sigma array. The cross sections are not computed for uniform grid of '
                              'pressures and temperatures')
                exit()

            sigma_3d = np.zeros((len(templist), len(preslist), self.nspecgrid))
            for idxtemp, valtemp in enumerate(templist):
                for idxpres, valpres in enumerate(preslist):
                    sigma_3d[idxtemp,idxpres,:] = self.readABSfiles(extfilelist=[absfilelist[valtemp][valpres]], # load only one file
                                                                    extpath=molpath,
                                                                    interpolate2grid=True,
                                                                    outputwavegrid = False)[0]

            # todo interpolation happens during fitting
            #build new temperature and pressure grids
            # tempgrid = np.arange(np.min(templist), np.max(templist), tempstep).tolist()
            # tempgrid.append(np.max(templist))
            # presgrid = np.arange(np.min(preslist), np.max(preslist), presstep).tolist()
            # presgrid.append(np.max(preslist))
            #
            # sigma_3d_reinterp = np.zeros((len(tempgrid), len(presgrid), self.nspecgrid))
            # for i in range(self.nspecgrid): # loop  wavelengths
            #     sigmainterp = interpolate.interp2d(preslist, templist, sigma_3d[:,:,i], kind='linear')
            #     sigma_3d_reinterp[:,:,i] = sigmainterp(tempgrid, presgrid).transpose()
            # moldict[molecule] = sigma_3d_reinterp

            moldict[molecule] = sigma_3d

        # sort lists
        templist = np.sort(templist)
        preslist = np.sort(preslist)

        # build sigma dictionary (in temperature and pressure)
        sigma_dict = {}
        for idxtemp, valtemp in enumerate(templist):
            if not valtemp in sigma_dict:
                sigma_dict[valtemp] = {}
            for idxpres, valpres in enumerate(preslist):
                if not valpres in sigma_dict[valtemp]:
                    sigma_dict[valtemp][valpres] = np.zeros((len(mollist), self.nspecgrid))
                j = 0
                for molecule in mollist:
                    sigma_dict[valtemp][valpres][j,:] = moldict[molecule][idxtemp,idxpres,:]
                    j += 1

        return sigma_dict, templist, preslist

    #@profile
    def readABSfiles(self,extfilelist=None, extpath= None,interpolate2grid=True,outputwavegrid=False):
    #reading in all absorption coefficient files and interpolating them to wavelength grid
    #
    # if filename = None, the absorption coefficient files will be read from the parameter file
    # if filename = [ascii list], absorption coefficient files will be read from user specified list
    # if outputwavegrid = True, it takes the first column of the ABS file and outputs as separate wavelength grid


        if extfilelist is None:
            #reading in list of abs files from parameter file
            if isinstance(self.params.in_abs_files,str):
                absfiles = genfromtxt(StringIO(self.params.in_abs_files),delimiter=",",dtype="S")
                try:
                    abslist = [absfiles.item()] #necessary for 0d numpy arrays
                except ValueError:
                    abslist = absfiles
            else:
                abslist = self.params.in_abs_files

            #checking if number of gasses in parameters file is consistent with gass columns in atm file
            if self.params.in_use_ATMfile and len(abslist) != self.ngas:
                # print len(abslist)
                # print self.ngas
                raise IOError('Number of gasses in .atm file incompatible with number of .abs files specified in parameters file')
                exit()

            out, wave = self.__readABSfiles_sub(path=self.params.in_abs_path, filelist=abslist,
                                                interpolate2grid=interpolate2grid, num=self.ngas)

        else:
            out, wave = self.__readABSfiles_sub(path=extpath, filelist=extfilelist,
                                                interpolate2grid=interpolate2grid,num=len(extfilelist))

        if outputwavegrid:
            return out, wave
        else:
            return out


    #@profile
    def __readABSfiles_sub(self, path, filelist, interpolate2grid, num):

        if interpolate2grid:
            out = np.zeros((num,self.nspecgrid))
            wave = np.transpose(self.readfile(os.path.join(path, filelist[0]), interpolate=True)[:,0])
        else:
            tmp = self.readfile(os.path.join(path, filelist[0]), interpolate=False)
            ABSsize = len(tmp[:,0])
            out = np.zeros((num,ABSsize))
            wave = np.transpose(tmp[:,0])

        for i in range(num):
            # print filelist[i]
            out[i,:] = np.transpose(self.readfile(os.path.join(path, filelist[i]), interpolate=interpolate2grid)[:,1])* 1e-4 #converting cm^2 to m^2

        return out, wave

    #@profile
    def set_ABSfile(self,path=None,filelist=None,interpolate = False,temperature=None):
        # @todo This should appear in the Preselector! CHECK!

        #manually overwrites absorption coefficients from new file
        #input path needs to be given and list of filename strings
        if path is None:
            extpath = self.params.in_abs_path
        if filelist is None:
            raise IOError('No input ABS file specified')
        if temperature is None:
            temperature = self.params.planet_temp
        sigma_array,self.specgrid = self.readABSfiles(extpath=path,
                                                      extfilelist=filelist,
                                                      interpolate2grid=interpolate,
                                                      outputwavegrid=True)
        self.sigma_dict = {}
        self.sigma_dict['tempgrid'] = [int(temperature)]
        self.sigma_dict[int(temperature)] = sigma_array
        self.nspecgrid = len(self.specgrid)
        
        
    #@profile
    def readATMfile(self):
    #reads in .atm file
        try:
            out = np.loadtxt(self.params.in_atm_file)
        except ValueError:
            out = np.loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
#         OUT[:,2] *= 1000. #converting from km to m

        out = out[np.argsort(out[:,2]),:]

        return out[:,0:3],np.transpose(out[:,3:])

    #@profile
    def build_rayleigh_sigma(self):

        # write rayleigh scattering cross sections to files

        logging.info('Write rayleigh scattering cross sections to files')

        gases = self.params.planet_molec + self.params.planet_inactive_gases

        for gasname in gases:

        # get the refractive index of a given gas gasname at a specific wavelength wl
        # Formulae taken from Allen Astrophysical Quantities if not otherwise specified

            n_formula = True # assume we have a formula for the refractive index of gasname
            king = 1 # King correction factor
            ns = 0   # refractive index
            wl = np.linspace(0.3, 30, 1000) # wavelengths in micron

            if gasname == 'He':
                ns = 1 + 0.01470091/(423.98-wl**-2) # C. R. Mansfield and E. R. Peck. Dispersion of helium, J. Opt. Soc. Am. 59, 199-203 (1969)
                # king is one for He
            elif gasname == 'H2':
                ns = 1 + 13.58e-5 * (1. + 7.52e-3 / wl**2)
                delta = 0.035 # from Morgan old code..
                king = (6.+3.*delta)/(6.-7.*delta) # Bates (1984)
            elif gasname == 'N2':
                ns = 1. + (6498.2 + (307.43305e12)/(14.4e9 - (1/(wl*1.e-4))**2))/1.e8
                #ns = 1 + 29.06e-5 * (1. + 7.7e-3 / wl**2)
                king = 1.034  #  Bates (1984)
            elif gasname == 'O2':
                ns = 1 + 1.181494e-4 + 9.708931e-3/(75.4-wl**-2) #  J. Zhang, Z. H. Lu, and L. J. Wang. Appl. Opt. 47, 3143-3151 (2008)
                king = 1.096
            elif gasname == '12C-16O2':
                #A. Bideau-Mehu, Y. Guern, R. Abjean and A. Johannin-Gilles. Interferometric determination of the refractive index of carbon dioxide in the ultraviolet region, Opt. Commun. 9, 432-434 (1973)
                ns = 1 + 6.991e-2/(166.175-wl**-2)+1.44720e-3/(79.609-wl**-2)+6.42941e-5/(56.3064-wl**-2)+5.21306e-5/(46.0196-wl**-2)+1.46847e-6/(0.0584738-wl**-2)
                king = 1.1364 #  Sneep & Ubachs 2005
            elif gasname == '12C-1H4':
                ns = 1 + 1.e-8*(46662. + 4.02*1.e-6*(1/(wl*1.e-4))**2)
            elif gasname == '12C-16O':
                ns = 1 + 32.7e-5 * (1. + 8.1e-3 / wl**2)
                king = 1.016  #  Sneep & Ubachs 2005
            elif gasname == '14N-1H3':
                ns = 1 + 37.0e-5 * (1. + 12.0e-3 / wl**2)
            elif gasname == '1H2-16O':
                ns_air = (1 + (0.05792105/(238.0185 - wl**-2) + 0.00167917/(57.362-wl**-2))) # P. E. Ciddor. Appl. Optics 35, 1566-1573 (1996)
                ns = 0.85 * (ns_air - 1.) + 1  # ns = 0.85 r(air) (Edlen 1966)
                delta = 0.17 # Marshall & Smith 1990
                king = (6.+3.*delta)/(6.-7.*delta)
            else:
                # this sets sigma_R to zero for all other gases
                n_formula = False
                logging.warning('There is no formula for the refractive index of %s. Cannot compute the cross section' % gasname)

            if n_formula: # only if the refractive index was computed
                Ns = 2.6867805e25 # in m^-3
                sigma_R =  ( 24.*np.pi**3/ (Ns**2) ) * (((ns**2 - 1.)/(ns**2 + 2.))**2) * king / (wl * 1.e-6)**4
                filename = os.path.join('Input', 'rayleigh', '%s.rayleigh' % gasname)
                logging.info('Saving cross section for %s to %s' % (gasname, filename))
                np.savetxt(filename, np.vstack((wl, sigma_R)).transpose())

    #@profile
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

    #@profile
    def get_molecular_weight(self, gasname):

        if gasname == 'He':
            mu = 4.
        elif gasname == 'H2':
            mu = 2.
        elif gasname == 'N2':
            mu = 28.
        elif gasname == 'O2':
            mu = 32.
        elif gasname == '12C-16O2':
            mu = 44.
        elif gasname == '12C-1H4':
            mu = 16.
        elif gasname == '12C-16O':
            mu = 28.
        elif gasname == '14N-1H3':
            mu = 17.
        elif gasname == '1H2-16O':
            mu = 18.
        elif gasname == 'C2H2':
            mu = 26.04
        elif gasname == '1H2-32S':
            mu = 34.0809
        else:
            mu = 0

        return mu * AMU