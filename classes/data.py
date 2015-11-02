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

np.set_printoptions(threshold=numpy.nan)


#loading taurex license manager. Only loaded in data class
import license
from license import *

try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
except:
    MPIrank = 0
    pass

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

        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # set observed spectrum specific variables (only if spectrum is provided)
            self.obs_wlgrid = self.obs_spectrum[:,0] # wavegrid in micron
            self.obs_nwlgrid = len(self.obs_spectrum[:,0]) # number of datapoints in spectrum
            self.obs_binwidths = self.obs_spectrum[:,3]   if shape(self.obs_spectrum)[1] == 4 else None # bin widths

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
            if self.obs_binwdiths == None:
                # if bin widths are *not* provided in the input spectrum
                bin_up -=  (self.obs_wlgrid[-1]-self.obs_wlgrid[-2])/2.
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

        # create wavenumber grid of internal model using numin, numax and delta wavenumber provided in param file
        # NOTE THAT DELTA NU MUST CORRESPOND TO THE CORRECT DELTA NU IN THE INPUT CROSSECTIONS
        # todo: this will change if we use non uniform grids. An external file with the grid should be provided...
        self.int_wngrid = np.arange(numin, numax, self.params.in_abs_dnu)
        self.int_dwngrid = np.diff(self.int_wngrid) # not very useful... this is an array of equal dnu
        self.int_nwngrid = len(self.int_wngrid) # number of points in the grid

        if isinstance(self.obs_spectrum, (np.ndarray, np.generic)):
            # calculate spectral binning grid
            self.obs_wngrid = 10000./self.obs_wlgrid
            self.intsp_bingrid, self.intsp_bingrididx = self.get_specbingrid(self.obs_wlgrid, self.int_wngrid, self.obs_binwidths)

        #reading in atmospheric profile file
        if self.params.in_use_ATMfile:
            logging.info('Using .atm file to generate atmospheric TP profile')
            self.pta, self.X = self.readATMfile() # pta = pressure, temp, alt; X = mixing ratios of molecules
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])
        else:
            self.nlayers = self.params.tp_atm_levels
            self.ngas = len(self.params.planet_molec)

        if MPIrank == 0:
            # compute rayleigh scattering cross sections and save them to a file
            self.build_rayleigh_sigma()
        self.sigma_R = self.get_rayleigh_sigma() # load Rayleigh scattering cross sections

        # Absorption cross sections
        if self.params.gen_run_gui:
            # if running in GUI mode, the dictionary is loaded from a file
            # todo needs some testing now...
            import pickle
            filename = 'Input/sigma_dict.pickle'
            if os.path.isfile(filename):
                # sigma dict has been presaved with the input parameters (lambda = 0.4-20, resolution = 500)
                logging.info('Loading sigma dictionary from file')
                self.sigma_dict = pickle.load(open(filename, 'rb'))
            else:
                logging.info('Creating sigma file')
                self.sigma_dict  = self.build_sigma_dic(tempstep=params.in_tempres)

        else:
            # reading in absorption coefficient data
            if self.params.in_use_P_broadening:
                # todo Pressure broadening: work in progress
                # self.sigma_dict_pres, self.sigma_templist, self.sigma_preslist  = self.build_sigma_dic_pressure()
                # tmpdb = self.sigma_dict_pres, self.sigma_templist, self.sigma_preslist
                # import pickle
                # pickle.dump(tmpdb, open('sigma_dict_press2.db', 'wb'))
                import pickle
                self.sigma_dict_pres, self.sigma_templist, self.sigma_preslist = pickle.load(open('sigma_dict_press2.db'))
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
    def get_specbingrid(self, wavegrid, specgrid, binwidths=None):
        #function calculating the bin boundaries for the data 
        #this is used to bin the internal spectrum to the data in fitting module

        if not isinstance(binwidths, (np.ndarray, np.generic)):
            bingrid =[]
            bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
            for i in range(len(wavegrid)-1):
                bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
            bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
            bingrid_idx = numpy.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
        else:
            # this bingrid is actually useless, as it doesn't allow for gaps in the data
            bingrid = []
            for i in range(len(wavegrid)):
                bingrid.append(wavegrid[i]-binwidths[i]/2.)
                bingrid.append(wavegrid[i]+binwidths[i]/2.)

            # build bin grid index array (an index for each model datapoint)
            bingrid_idx = np.empty(len(specgrid))
            bingrid_idx[:] = np.NaN
            for i in range(len(specgrid)):
                for j in range(len(wavegrid)):
                    if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                        bingrid_idx[i] = j+1
                        break

        return bingrid, bingrid_idx

    def get_specgrid_from_crosssections(self, lambda_min=0.1, lambda_max=20.0):
        # return the wavelength grid of the absorbion cross sections
        # assume that the cross sections have a uniform wl grid

        if self.params.in_use_P_broadening:
            path = self.params.in_abs_path_P
        else:
            path = self.params.in_abs_path

        files = glob.glob(os.path.join(path, '*.txt'))
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
            SED = libem.black_body(self.int_wngrid,self.params.star_temp) #@todo bug here? not multiplied by size of star 4piRs^2
#             SED *= self.params.star_radius**2 * np.pi * 4.

#             SED *= self.params.star_radius**2 * np.pi * 4.
        else:
            # finding closest match to stellar temperature in parameter file
            [tmpselect, idx] = libgen.find_nearest(tmpind, self.params.star_temp)
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
#         SED = libem.black_body(self.int_wngrid,self.params.star_temp)
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

        if interpolate:
            # inteprolate to wavenumber grid
            interpflux = interp(self.int_wngrid, out[:,0], out[:,1])
            out = transpose(vstack((self.int_wngrid, interpflux)))

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

            sigtmp = np.zeros((len(absfiles), self.int_nwngrid))

            for i in range(len(absfiles)): # loop through xsecs for given molecule


                # read cross section file
                xsec = np.loadtxt(os.path.join(self.params.in_abs_path, molecule, absfiles[i]))

                # paste the input cross section absfiles[i] into the wavenumber grid of the sigma array

                # the following takes care of stripping out the bits absfiles[i] that fall outside the internal
                # wavenumber grid self.int_wngrid. Also, for regions of self.int_wngrid that are not covered by
                # absfiles[i] we assume sigma = 0
                idmin_xsec = 0
                idmin_sigdic = 0
                if np.min(xsec[:,0]) < np.min(self.int_wngrid):
                    idmin_xsec = np.where(np.abs(xsec[:,0] - np.min(self.int_wngrid)) < 1e-6)[0][0]
                else:
                    idmin_sigdic = np.where(np.abs(self.int_wngrid - np.min(xsec[:,0])) < 1e-6)[0][0]

                idmax_xsec = len(xsec[:,0]) - 1
                idmax_sigdic = self.int_nwngrid - 1
                if np.max(xsec[:,0]) > np.max(self.int_wngrid):


                    idmax_xsec = np.where(np.abs(xsec[:,0] - np.max(self.int_wngrid)) < 1e-6)[0][0]
                else:
                    idmax_sigdic = np.where(np.abs(self.int_wngrid - np.max(xsec[:,0])) < 1e-6)[0][0]
                #
                #print idmin_sigdic, idmax_sigdic, idmin_xsec, idmax_xsec
                # print shape(xsec[:,0])
                # print np.min(xsec[:,0]), np.max(xsec[:,0])
                # print shape(sigtmp[i,idmin_sigdic:idmax_sigdic]), shape(xsec[:,1][idmin_xsec:idmax_xsec])


                sigtmp[i,idmin_sigdic:idmax_sigdic] =  xsec[:,1][idmin_xsec:idmax_xsec] * 1e-4 # converting cm^2 to m^2
            sigshape = np.shape(sigtmp)
            moldict[molecule]['sigma'] = sigtmp
            moldict[molecule]['templist'] = templist
            moldict[molecule]['tempmin'] = tempmin
            moldict[molecule]['tempmax'] = tempmax

        # setting up new temperature grid
        temp_totmin = np.max(tempmin)
        temp_totmax = np.min(tempmax)
        tempgrid = [i for i in range(int(temp_totmin),int(temp_totmax),int(tempstep))]

        tempgrid.append(int(temp_totmax))
        moldict['tempgrid'] = tempgrid

        # interpolating absorption cross sections to new temperature grid
        for molecule in mollist:
            interpsigma = np.zeros((len(tempgrid),sigshape[1]))
            for i in range(sigshape[1]):
                # todo: should rather use the interpolation proposed by Hill et al. Alternatively, we can interpolate in the forward model, in the C code... Need to check what's faster (array slicing vs live interpolation).
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
            molpath = os.path.join(self.params.in_abs_path_P)
            absfilelist, templist_tmp, preslist_tmp = libgen.find_absfiles_pressure(molpath, molecule)

            if (templist == None or templist == templist_tmp) and (preslist == None or preslist == preslist_tmp):
                templist = templist_tmp
                preslist = preslist_tmp
            else:
                logging.error('Cannot build sigma array. The cross sections are not computed for uniform grid of '
                              'pressures and temperatures')
                exit()

            sigma_3d = np.zeros((len(templist), len(preslist), self.int_nwngrid-1))
            for idxtemp, valtemp in enumerate(templist):
                for idxpres, valpres in enumerate(preslist):
                    xsec =  self.readfile(os.path.join(molpath, molecule, absfilelist[valtemp][valpres]))
                    sigma_3d[idxtemp,idxpres,:] =  xsec[:,1][np.logical_and(xsec[:,0]>np.min(self.int_wngrid), xsec[:,0]<(np.max(self.int_wngrid)))] * 1e-4 #converting cm^2 to m^2

            # todo interpolation happens during fitting
            #build new temperature and pressure grids
            # tempgrid = np.arange(np.min(templist), np.max(templist), tempstep).tolist()
            # tempgrid.append(np.max(templist))
            # presgrid = np.arange(np.min(preslist), np.max(preslist), presstep).tolist()
            # presgrid.append(np.max(preslist))
            #
            # sigma_3d_reinterp = np.zeros((len(tempgrid), len(presgrid), self.int_nwngrid))
            # for i in range(self.int_nwngrid): # loop  wavelengths
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
                    sigma_dict[valtemp][valpres] = np.zeros((len(mollist), self.int_nwngrid-1))
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
            out = np.zeros((num,self.int_nwngrid))
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
        sigma_array,self.int_wngrid = self.readABSfiles(extpath=path,
                                                      extfilelist=filelist,
                                                      interpolate2grid=interpolate,
                                                      outputwavegrid=True)
        self.sigma_dict = {}
        self.sigma_dict['tempgrid'] = [int(temperature)]
        self.sigma_dict[int(temperature)] = sigma_array
        self.int_nwngrid = len(self.int_wngrid)


    #@profile
    def readATMfile(self):
    #reads in .atm file
        try:
            out = np.loadtxt(self.params.in_atm_file)
        except ValueError:
            out = np.loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
        out[:,2] *= 1000. #converting from km to m

        out = out[np.argsort(out[:,2]),:]

        return out[:,0:3],np.transpose(out[:,3:])

    #@profile
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