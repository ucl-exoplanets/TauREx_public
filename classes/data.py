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

class data(base):

    # @todo general comment: maybe move some stuff from this class to the tp_profile, if appropriate

    def __init__(self, params, spectrum=None):

        logging.info('Initialising data object')

        #checking taurex license
        license_manager().run()
        
        self.params = params
        self.KBOLTZ = 1.380648813e-23
        self.G      = 6.67384e-11
        
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

        #setting planetary surface gravity
        self.planet_grav = self.get_surface_gravity()
 
        
        #calculating atmospheric scale height
        self.scaleheight = self.get_scaleheight() # @todo needed? Check. Scaleheight is more specific to tp_profile class

        #reading in atmospheric profile file
        #or generating it from parameter file values
        #pta = pressure, temp, alt
        #X = mixing ratios of molecules
        if self.params.in_use_ATMfile: # @todo maybe move to tp_profile?
            self.pta, self.X = self.readATMfile()
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])
        else:
            #if no pta file provided, pta will be calculated by the profile class
            self.nlayers = self.params.tp_atm_levels
            self.ngas    = size(self.params.planet_molec)
#             self.pta     = self.setup_pta_grid() #if not read from file, pta will be calculated by profile class
            self.X       = self.set_mixing_ratios() #mixing ratios are read from parameter file or preselector class @todo move to atmosphere class?


        #setting up dictionary with atmosphere parameters
        self.atmosphere = self.init_atmosphere() # @todo this needs to be improved --> move to tp_profile
        
        #reading in absorption coefficient data 
#         self.sigma_array = self.readABSfiles()
        self.sigma_dict = self.build_sigma_dic(tempstep=params.in_tempres) # @todo move to tp_profile?

        #reading in other files if specified in parameter file
        if params.in_include_rad:
            self.rad = self.readfile(self.params.in_rad_file, interpolate=True)
        if params.in_include_cia:
            self.cia = self.readfile(self.params.in_cia_file, interpolate=True)
#         if params.in_include_cld:
#             self.cld = self.readfile(self.params.in_cld_file,interpolate=True)
            
        #reading in Phoenix stellar model library (if emission is calculated only)
        if self.params.gen_type == 'emission' or self.params.fit_emission:
            self.F_star = self.get_star_SED() #@todo there is most certainly a bug there. think units are ergs/s/cm^2 at the moment 

    #class functions

    def init_atmosphere(self, mu=0.0, def_mu=2.3):
        #initialising atmosphere dictionary
        ATM = {}
        ATM['mol'] = {}
        ATM['info'] = {}
        ATM['info']['def_mu']  = def_mu #default atmos 85% H2, 15% He --> mu~2.3
        ATM['info']['mu']      = mu
        ATM['info']['nmol']    = 0      #number of molecule index
        
        return ATM

    def get_surface_gravity(self):
        #calculate surface gravity of planet using Rp and Mp
        
        return (self.G * self.params.planet_mass) / (self.params.planet_radius**2)

    def get_scaleheight(self,T_aver=None,surf_g=None,mmw=None):
        #compute scaleheight of atmosphere
        if T_aver is None:
            T_aver = self.params.planet_temp
        if surf_g is None:
            surf_g = self.planet_grav
        if mmw is None:
            mmw = self.params.planet_mu
             
 
        return (self.KBOLTZ*T_aver)/(mmw*surf_g)

        

    def get_specgrid(self,R=5000,lambda_min=0.1,lambda_max=20.0):
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

        # figure()
        # plot(specgrid)
        # figure()
        # plot(delta_lambda)
        # show()
        # exit()

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

    def set_mixing_ratios(self):

        #setting up mixing ratio array from parameter file inputs

        mixing = self.params.planet_mixing

        #checking if number of mixing ratios = number of gasses
        X = zeros((self.ngas,self.nlayers))
        
        if self.params.pre_run:
            X += 1e-4
        else:
            if size(mixing) is not self.ngas:
                raise IOError('Wrong  number of mixing ratios to molecules in parameter file')
                exit()

            # X = np.tile(mixing, [self.nlayers, 1]).transpose()
            for i in range(self.ngas):
                X[i,:] += float(mixing[i])
        return X
        
        

    def add_molecule(self,NAME, WT, RAD, RIDX,FRAC):
    #adding a molecule to the atmosphere dictionary
        self.atmosphere['mol'][NAME] = {}
        self.atmosphere['mol'][NAME]['weight']  = WT   #relative molecular weight (amu)
        self.atmosphere['mol'][NAME]['radius']  = RAD  #molecular radius (m)
        self.atmosphere['mol'][NAME]['ridx']    = RIDX #refractive index
        self.atmosphere['mol'][NAME]['frac']    = FRAC #fraction of total composition
        self.atmosphere['info']['nmol']         += 1   #increasing molecule count
        #updating mean molecular weight of atmosphere
        self.atmosphere['info']['mu'] = self.get_mean_molweight()
        
    def get_mean_molweight(self):
    #returning mean molecular weight of atmosphere
        MMW = 0
        for m,M in enumerate(self.atmosphere['mol']):
            MMW += (self.atmosphere['mol'][M]['weight'] * self.atmosphere['mol'][M]['frac'])
        
        return MMW


    def set_ABSfile(self,path=None,filelist=None,interpolate = False,temperature=None):
    #manually overwrites absorption coefficients from new file
    #input path needs to be given and list of filename strings
        if path is None:
            extpath = self.params.in_abs_path
        if filelist is None:
            raise IOError('No input ABS file specified')
        if temperature is None:
            temperature = self.params.planet_temp
        
        sigma_array,self.specgrid = self.readABSfiles(extpath=path,extfilelist=filelist, interpolate2grid=interpolate,outputwavegrid=True)
        self.sigma_dict = {}
        self.sigma_dict['tempgrid'] = [int(temperature)]
        self.sigma_dict[int(temperature)] = sigma_array
        self.nspecgrid = len(self.specgrid)



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
                print 'WARNING: Stellar temp. in .par file exceeds range ', min(tmpind), '-', max(tmpind), 'K'
                print 'Using black-body approximation instead'
            self.star_blackbody = True
            SED = libem.black_body(self.specgrid,self.params.star_temp) #@todo bug here! not multiplied by size of star 4piRs^2
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
            SED_raw *= 10.
            SED = np.interp(self.specgrid, SED_raw[:,0], SED_raw[:,1])
    
        return SED
    


    def build_sigma_dic(self,tempstep=50):
        #building temperature dependent sigma_array
        
#        try:
#            mollist = [self.params.planet_molec.item()] #necessary for 0d numpy arrays
#        except ValueError:
#            mollist = self.params.planet_molec

        mollist = self.params.planet_molec

#         print self.params.planet_molec
        
        tempmin = []
        tempmax = []
                
        moldict = {}
        
        
        for molecule in mollist:

            logging.info('Load sigma array for molecule %s' % molecule)

            moldict[molecule] ={}
            
            absfiles, templist = libgen.find_absfiles(self.params.in_abs_path,molecule)
            tempmax.append(np.max(templist))
            tempmin.append(np.min(templist))
            
            sigtmp = self.readABSfiles(extfilelist=absfiles, extpath=self.params.in_abs_path, interpolate2grid=True, outputwavegrid = False)
            
            sigshape = np.shape(sigtmp)
            
            moldict[molecule]['sigma'] = sigtmp
            moldict[molecule]['templist'] = templist
            moldict[molecule]['tempmin'] = tempmin
            moldict[molecule]['tempmax'] = tempmax
            
            
#             figure()
#             plot(np.transpose(sigtmp))
#             show()

        #setting up new temperature grid
        temp_totmin = np.max(tempmin)
        temp_totmax = np.min(tempmax)
        tempgrid = [i for i in range(int(temp_totmin),int(temp_totmax),int(tempstep))]
        tempgrid.append(int(temp_totmax))
        moldict['tempgrid'] = tempgrid 
        
#         print temp_totmax, temp_totmin, tempgrid
        
        #interpolating absorption cross sections to new temperature grid
        for molecule in mollist:
            interpsigma = zeros((len(tempgrid),sigshape[1]))
            
            for i in range(sigshape[1]):
                interpsigma[:,i] = np.interp(tempgrid,moldict[molecule]['templist'],moldict[molecule]['sigma'][:,i]) #linear interpolation. to be changed with hill et al method 
            
            moldict[molecule]['interpsigma'] = interpsigma
            
#             print np.shape(interpsigma)
#             figure()
#             plot(np.transpose(interpsigma))

        #building final sigma_array dictionary
        sigma_dict = {}
        sigma_dict['tempgrid'] = tempgrid
        for i in range(len(tempgrid)):
            sigma_dict[tempgrid[i]] = {}
            
            sigma_array = zeros((len(mollist),sigshape[1]))
            j = 0
            for molecule in mollist:
                sigma_array[j,:] = moldict[molecule]['interpsigma'][i,:]
                j += 1
            sigma_dict[tempgrid[i]] = sigma_array
            
#             print tempgrid[i], type(tempgrid[i])
#             figure()
#             plot(np.transpose(sigma_dict[tempgrid[i]]))
#             show() 
            
        
#         figure()
#         plot(np.transpose(sigma_dict[600][0,:]),'b')
#         plot(np.transpose(sigma_dict[1600][0,:]),'r')
#         
#         show()

        return sigma_dict



    def readATMfile(self):
    #reads in .atm file
        try:
            OUT = loadtxt(self.params.in_atm_file)
        except ValueError:
            OUT = loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
#         OUT[:,2] *= 1000. #converting from km to m
        
        OUT = OUT[argsort(OUT[:,2]),:]
        
        return OUT[:,0:3],transpose(OUT[:,3:])

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

            OUT, WAVE = self.__readABSfiles_sub(path=self.params.in_abs_path,filelist=abslist, interpolate2grid=interpolate2grid,num=self.ngas)

        else:
            OUT, WAVE = self.__readABSfiles_sub(path=extpath,filelist=extfilelist, interpolate2grid=interpolate2grid,num=len(extfilelist))

        if outputwavegrid:
            return OUT, WAVE
        else:
            return OUT
        

    def __readABSfiles_sub(self,path, filelist, interpolate2grid,num):
        if interpolate2grid:
            OUT = zeros((num,self.nspecgrid))
            WAVE = transpose(self.readfile(path+filelist[0], interpolate=True)[:,0])
        else:
            tmp = self.readfile(path+filelist[0], interpolate=False)
            ABSsize = len(tmp[:,0])
            OUT = zeros((num,ABSsize))
            WAVE = transpose(tmp[:,0])

        for i in range(num):
            # print filelist[i]
            OUT[i,:] = transpose(self.readfile(path+filelist[i], interpolate=interpolate2grid)[:,1])* 1e-4 #converting cm^2 to m^2

        return OUT, WAVE


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

#         figure()
#         plot(out[:,0], out[:,1])

        #interpolating to specgrid
        if interpolate:
#             interpflux = interp1d(out[:,0],out[:,1],axis=0,kind='cubic')(self.wavegrid)
            interpflux = interp(self.specgrid, out[:,0], out[:,1])
#             print interpflux
            out = transpose(vstack((self.specgrid, interpflux)))
        # print 'ble',np.shape(out)
#         plot(out[:,0], out[:,1], c='r')
#         show()
        return out

        