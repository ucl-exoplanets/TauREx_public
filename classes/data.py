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
import library_general as libgen
import library_emission as libem
import logging

#loading taurex license manager. Only loaded in data class
import license
from license import *

# some constants
KBOLTZ = 1.380648813e-23
G = 6.67384e-11

class data(base):

    def __init__(self, params, spectrum=None):

        logging.info('Initialising data object')

        #checking taurex license
        license_manager().run()
        
        self.params = params

        #converting absorption cross-sectinos from cm^-1 to microns
        if params.in_convert2microns:
            libgen.convert2microns(params.in_abs_path+'*')
        
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
            self.specgrid, self.dlamb_grid = self.get_specgrid(R=self.params.gen_spec_res,
                                                               lambda_min=self.params.gen_wavemin,
                                                               lambda_max=self.params.gen_wavemax)
        else:
            # user wavelength range from observed spectrum, but increase model wavelength range by one data spectral bin
            # in blue and red. This ensures correct model binning at the edges
            bin_up = self.wavegrid[-1]-self.wavegrid[-2]
            bin_low = self.wavegrid[1]-self.wavegrid[0]
            self.specgrid, self.dlamb_grid = self.get_specgrid(R=self.params.gen_spec_res,
                                                               lambda_min=self.wavegrid[0]-bin_low,
                                                               lambda_max=self.wavegrid[-1]+bin_up)
        self.nspecgrid = len(self.specgrid)

        if isinstance(self.spectrum, (np.ndarray, np.generic)):
            #calculating spectral binning grid, only if observed spectrum is provided
            self.spec_bin_grid, self.spec_bin_grid_idx = self.get_specbingrid(self.wavegrid)
            self.n_spec_bin_grid= len(self.spec_bin_grid)

        #calculating atmospheric scale height
        #self.scaleheight = self.get_scaleheight() # @todo needed? Check. Scaleheight is more specific to tp_profile class

        #reading in atmospheric profile file
        #or generating it from parameter file values
        #pta = pressure, temp, alt
        #X = mixing ratios of molecules
        if self.params.in_use_ATMfile:
            self.pta, self.X = self.readATMfile()
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])

        #reading in absorption coefficient data
        self.sigma_dict  = self.build_sigma_dic(tempstep=params.in_tempres)

        # #reading in other files if specified in parameter file
        if params.in_include_rad:
            self.rad = self.readfile(self.params.in_rad_file, interpolate=True)
        if params.in_include_cia:
            self.cia = self.readfile(self.params.in_cia_file, interpolate=True)
        if params.in_include_cld:
            self.cld = self.readfile(self.params.in_cld_file, interpolate=True)


        #reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.F_star = self.get_star_SED() #@todo there is most certainly a bug there. think units are ergs/s/cm^2 at the moment 

        # list of all molecules for which we have cross sections
        self.all_absorbing_gases = ['1H2-16O', '1H-12C-14N', '12C-1H4', '12C-16O2', '12C-16O', '14N-1H3',
                                    '28Si-16O', '48Ti-16O', '51V-16O']

        # list of all inactive gases we take care of
        self.all_inactive_gases = ['He', 'H2', 'N2']

        # dictionary with parameters (e.g. mean molecular weight, refractive index, etc) for each gas (absorbers and inactive)
        self.all_gases_properties = get_gases_properties()

    #class functions

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
    
    def get_specbingrid(self,wavegrid):
        #function calculating the bin boundaries for the data 
        #this is used to bin the internal spectrum to the data in fitting module
        
        bingrid =[]
        bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
        for i in range(len(wavegrid)-1):
            bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
        bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
        
        bingrid_idx = numpy.digitize(self.specgrid,bingrid) #getting the specgrid indexes for bins
        

        return bingrid, bingrid_idx 



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
                logging.warning('WARNING: Stellar temp. in .par file exceeds range %.1f - %.1f K. '
                                'Using black-body approximation instead' % (min(tmpind), max(tmpind)))
            self.star_blackbody = True
            SED = libem.black_body(self.specgrid,self.params.star_temp) #@todo bug here! not multiplied by size of star 4piRs^2
#             SED *= self.params.star_radius**2 * np.pi * 4.

#             SED *= self.params.star_radius**2 * np.pi * 4.
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = libgen.find_nearest(tmpind, self.params.star_temp)
            self.star_blackbody = False
            
            for file in fileindex: #this search is explicit due to compatibility issues with Mac and Linux sorting
                if np.int(file.split('/')[-1][3:8]) == np.int(tmpselect):
                    self.SED_filename = file

            #reading in correct file and interpolating it onto self.specgrid
            SED_raw = np.loadtxt(self.SED_filename, dtype='float', comments='#') #@todo bug here! not multiplied by size of star 4piRs^2 and units are wrong as well
            SED = np.interp(self.specgrid, SED_raw[:,0], SED_raw[:,1])
        return SED



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

        #interpolating to specgrid
        if interpolate:
            interpflux = interp(self.specgrid, out[:,0], out[:,1])
            out = transpose(vstack((self.specgrid, interpflux)))

        return out

    def build_sigma_dic(self,tempstep=50):

        #building temperature dependent sigma_array

        mollist = self.params.planet_molec

        tempmin = []
        tempmax = []

        moldict = {}

        for molecule in mollist:

            logging.info('Load sigma array for molecule %s' % molecule)

            moldict[molecule] ={}

            absfiles, templist = libgen.find_absfiles(self.params.in_abs_path,molecule)
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
                print len(abslist)
                print self.ngas
                raise IOError('Number of gasses in .atm file incompatible with number of .abs files specified in parameters file')
                exit()

            out, wave = self.__readABSfiles_sub(path=self.params.in_abs_path,filelist=abslist, interpolate2grid=interpolate2grid,num=self.ngas)

        else:
            out, wave = self.__readABSfiles_sub(path=extpath,filelist=extfilelist, interpolate2grid=interpolate2grid,num=len(extfilelist))

        if outputwavegrid:
            return out, wave
        else:
            return out


    def __readABSfiles_sub(self, path, filelist, interpolate2grid,num):
        if interpolate2grid:
            out = np.zeros((num,self.nspecgrid))
            wave = np.transpose(self.readfile(path+filelist[0], interpolate=True)[:,0])
        else:
            tmp = self.readfile(path+filelist[0], interpolate=False)
            ABSsize = len(tmp[:,0])
            out = np.zeros((num,ABSsize))
            wave = np.transpose(tmp[:,0])

        for i in range(num):
            # print filelist[i]
            out[i,:] = np.transpose(self.readfile(path+filelist[i], interpolate=interpolate2grid)[:,1])* 1e-4 #converting cm^2 to m^2

        return out, wave


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
        self.data.nspecgrid = len(self.specgrid)
        
        
    def readATMfile(self):
    #reads in .atm file
        try:
            out = np.loadtxt(self.params.in_atm_file)
        except ValueError:
            out = np.loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
#         OUT[:,2] *= 1000. #converting from km to m

        out = out[np.argsort(out[:,2]),:]

        return out[:,0:3],np.transpose(out[:,3:])

    def get_refractive_index(self, gasname, wl):

        # returns the refractive index of a given gas gasname at a specific wavelength wl
        # Formulae taken from Allen Astrophysical Quantities if not otherwise specified

        if gasname == 'He':
            r = 1 + 3.48e-5 * (1. + 2.3e-3 / wl**2)
        elif gasname == 'H2':
            r = 1 + 13.58e-5 * (1. + 7.52e-3 / wl**2)
        elif gasname == 'N2':
            r = 1 + 29.06-5 * (1. + 7.7e-3 / wl**2)
        elif gasname == 'O2':
            r = 1 + 26.63-5 * (1. + 5.07e-3 / wl**2)
        elif gasname == '12C-16O2':
            r = 1 + 43.9e-5 * (1. + 6.4e-3 / wl**2)
        elif gasname == '12C-1H4':
            r = 1.000441,
        elif gasname == '12C-16O':
            r = 1 + 32.7e-5 * (1. + 8.1e-3 / wl**2)
        elif gasname == '14N-1H3':
            r = 1 + 37.0e-5 * (1. + 12.0e-3 / wl**2)
        elif gasname == '1H2-16O':
            # r = 0.85 r(air) (Edlen 1966)
            # dispersion formula for air: P. E. Ciddor. Refractive index of air: new equations for the visible and near infrared, Appl. Optics 35, 1566-1573 (1996)
            r = 0.85 * (1 + (0.05792105/(238.0185 - wl**-2) + 0.00167917/(57.362-wl**-2)))
        else:
            r = 1

        return r

    def get_molecular_weight(self, gasname):

        if gasname == 'He':
            mw = 4.
        elif gasname == 'H2':
            mw = 2.
        elif gasname == 'N2':
            mw = 14.
        elif gasname == 'O2':
            mw = 32.
        elif gasname == '12C-16O2':
            mw = 44.
        elif gasname == '12C-1H4':
            mw = 16.
        elif gasname == '12C-16O':
            mw = 28.
        elif gasname == '14N-1H3':
            mw = 17.
        elif gasname == '1H2-16O':
            mw = 18.

        return mw