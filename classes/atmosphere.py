################################################
# class atmosphere (once called tp_profile)
#
# Reads in all relevant data and performs pre-processing
# such as sorting, grid-interpolations etc.
#
# Input: -data object
#        -parameter object (optional)
#
#
# Output: - T-P profiles
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013
#
################################################

#loading libraries
from base import base
import numpy as np
import pylab as pl
import logging

import library_general as libgen


# some constants
KBOLTZ = 1.380648813e-23
G = 6.67384e-11

class atmosphere(base):

    def __init__(self, data, params=None):

        logging.info('Initialising atmosphere object')

        # set some objects and variables
        if params is not None:
            self.params = params
        else:
            self.params = data.params
        self.data = data
        self.fit_transmission = self.params.fit_transmission
        self.fit_emission     = self.params.fit_emission

        #derived values @todo check planet_mu?
        self.scaleheight = self.get_scaleheight(self.params.planet_temp, data.planet_grav, self.params.planet_mu)

        # build PTA profile and mixing ratio array (pta = pressure, temp, alt; X = mixing ratios of molecules)
        if self.params.tp_var_atm:
            # atmospheric profile is fitted
            self.nlayers = int(self.params.tp_atm_levels)
            self.ngas = int(len(self.params.planet_molec))
            self.X = self.set_mixing_ratios() # mixing ratios are read from parameter file or set to 1e-4 if preselector = True
            self.pta = self.setup_pta_grid()

        elif self.params.in_use_ATMfile:
            # reading atmospheric profile from .atm file
            self.nlayers = len(self.pta[:,0])
            self.ngas = len(self.X[:,0])
            self.X = np.zeros((self.ngas,self.nlayers)) + (1e-5) #setting up initial mixing ratios
            self.pta, self.X = self.readATMfile()




        self.P = self.pta[:,0] # pressure array
        self.P_bar = self.P * 1.0e-5 #convert pressure from Pa to bar
        self.T = self.pta[:,1] # temperature array
        self.z = self.pta[:,2] # altitude array
        self.rho = self.get_rho(T=self.T, P=self.P)

        self.num_T_params = 3 #number of free temperature parameters @todo what is it ?

        if self.params.fit_emission or self.params.fit_transmission:
            self.bounds = self.get_prior_bounds()

        if self.params.fit_emission:
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(emission=True)

        if self.params.fit_transmission:
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(transmission=True)

        # reading in absorption coefficient data.
        self.sigma_dict = self.build_sigma_dic(tempstep=self.params.in_tempres)

        #reading in other files if specified in parameter file
        if self.params.in_include_rad:
            self.rad = self.data.readfile(self.params.in_rad_file, interpolate=True)
        if self.params.in_include_cia:
            self.cia = self.data.readfile(self.params.in_cia_file, interpolate=True)
#
#         self.fit_params[0] = 5e-6
#         self.fit_params[1] = 2e-7
#         self.fit_params[2] = 1400
#         self.fit_params[3] = 1400
#         self.fit_params[4] = 1200
#         self.fit_params[5] = 1e5
#         self.fit_params[6] = 100.0


#         self.bounds = [(0.0, 0.01), (0.0, 0.01), (1000.0, 1800.0), (1000.0, 1800.0), (1000.0, 1800.0), (50000.0, 500000.0), (50.0, 150.0)]
#         self.Tpriors = [(1000.0, 1800.0), (1000.0, 1800.0), (1000.0, 1800.0)]
#         self.Ppriors = [(50000.0, 500000.0), (50.0, 150.0)]
#         self.Xpriors = [(0.0, 0.01), (0.0, 0.01)]
#         print self.bounds
#
#
#         self.fit_params,self.fit_index, self.fit_count = self.setup_parameter_grid(fit_emission=True)
#         print self.fit_params
#         PARAMS2 = self.fit_params
#         PARAMS2[2] = 1400
#         PARAMS2[3] = 1400
#         PARAMS2[4] = 1200
# #
# #         print PARAMS2
# #
#         self.T,self.P,self.X = self.TP_profile(PARAMS=PARAMS2)
#         self.rho = self.get_rho(T=self.T,P=self.P)
#
#         pl.figure(104)
#         pl.plot(T,np.log(P))
#
#         rho = self.get_rho(T=T, P=P)
#         pl.figure(105)
#         pl.plot(rho)
#         pl.show()
#         exit()
#         pl.figure(200)
#         pl.plot(T,P)
#         pl.yscale('log')
#         pl.gca().invert_yaxis()
#         pl.show()
#         exit()

    #class methods

    def get_scaleheight(self,T,g,mu):
        return (KBOLTZ*T)/(mu*g)

    def setup_pta_grid(self, T=None):

        #calculate pressure, temperature, altitude grid

        max_p    = self.params.tp_max_pres
        n_scale  = self.params.tp_num_scale # thickness of atmosphere in number of atmospheric scale heights
        n_layers = self.nlayers

        # get new scale height if T is provided. Otherwise assume default (should be params.planet_temp)
        if T is not None:
            self.scaleheight = self.get_scaleheight(T[0],  self.data.planet_grav, self.params.planet_mu)

        max_z = n_scale * self.scaleheight

        #generatinng altitude-pressure array
        pta_arr = np.zeros((n_layers,3))
        pta_arr[:,2] = np.linspace(0,max_z,num=n_layers) # altitude
        pta_arr[:,0] = max_p * np.exp(-pta_arr[:,2]/self.scaleheight)
#         if T is not None:
#             PTA_arr[:,1] = T
#         else:
        pta_arr[:,1] = self.params.planet_temp

        return pta_arr

    def get_prior_bounds(self):
        #partially to be moved to parameter file i guess

        self.Xpriors = [self.params.fit_X_low, self.params.fit_X_up] #lower and upper bounds for column density
        self.Tpriors = [self.params.planet_temp - self.params.fit_T_low, self.params.planet_temp + self.params.fit_T_up] #lower and upper bounds for temperature

        #BE REALLY CAREFUL WITH THIS PARAMETER. THIS NEEDS SANITY CHECKING
        self.Ppriors = [[5e4,5e5],[50.0,150.0]] #lower and upper bounds for individual TP transistion (e.g. tropopause height, mesospehere height) in Pa

        #setting up bounds for downhill algorithm
        #this may be merged into setup_parameter_grid() later but unsure of how complex this is going to be right now
        bounds = []
        for i in xrange(self.ngas):
            bounds.append((self.Xpriors[0],self.Xpriors[1]))
        if self.fit_transmission:
            bounds.append((self.Tpriors[0],self.Tpriors[1]))
        if self.fit_emission:
            for i in xrange(self.num_T_params):
                bounds.append((self.Tpriors[0],self.Tpriors[1]))
            for i in xrange(self.num_T_params-1):
                bounds.append((self.Ppriors[i][0],self.Ppriors[i][1]))

        return bounds



    def setup_parameter_grid(self, transmission=False, emission=False):

        # fit_params    [abundances (ngas), temperatures (num_T_params), pressures (num_T_params-1)]
        # fit_count     [no. abundances, no. temperatures, no. pressures]
        # fit_index     index to fitparams for slicing

        fit_params = []
        fit_count = []
        fit_index = []

        # setting up mixing ratios for individual gases
        Xmean = np.mean(self.Xpriors)
        cgas = 0
        for i in xrange(self.ngas):
            fit_params.append(Xmean)
            cgas += 1

        fit_count.append(cgas)

        #setting up temperature parameters
        T_mean = np.mean(self.Tpriors)
        num_T_params = self.num_T_params

        ctemp = 0; cpres = 0

        if transmission:
            ctemp +=1
            fit_params.append(self.params.planet_temp)

        if emission:
            for i in xrange(num_T_params):
                fit_params.append(T_mean)
                ctemp += 1
            for i in xrange(num_T_params-1):
                fit_params.append(np.mean(self.Ppriors[i]))
                cpres += 1

        fit_count.append(ctemp)
        fit_count.append(cpres)

        cumidx = fit_count[0]
        fit_index.append(cumidx)

        for i in xrange(1,len(fit_count)):
            cumidx += fit_count[i]
            fit_index.append(cumidx)

        return fit_params, fit_index, fit_count


    def TP_profile(self, fit_params, T=None, P=None):

    # @todo do we need T and P ? maybe not

    # main function defining basic parameterised TP-profile from
    # PARAMS and INDEX. INDEX = [Chi, T, P]

        fit_index = self.fit_index
        fit_count = self.fit_count

        X_params = fit_params[:fit_index[0]]
        T_params = fit_params[fit_index[0]:fit_index[1]]
        P_params = fit_params[fit_index[1]:]

        #convert X into 2d arrays (previously done in fitting module but seems more appropriate here)
#         X   = np.zeros((self.ngas,self.nlayers))
        self.X[:] = 0.0
        for i in xrange(self.ngas):
            self.X[i,:] += X_params[i]

        # Recalculate PTA profile, based on new Temperature.
        if len(T_params) > 0:
            self.pta = self.setup_pta_grid(T_params)
            P = self.pta[:,0]
        elif P is None:
            P = self.P

        # probably not needed?
        if T is not None:
            return T, P, self.X

        # if we have more than 1 temperature in fit_params
        if fit_count[1] > 1:
            P_params =  [self.params.tp_max_pres] + list(P_params) + [np.min(P)]
            T_params = list(T_params) + [T_params[-1]]

            #creating linear T-P profile
            T = np.interp(np.log(P[::-1]), np.log(P_params[::-1]), T_params[::-1])
            return T[::-1], P, self.X

        if fit_count[1] == 1:
            T = T_params
            return T, P, self.X

    def get_rho(self, T=None, P=None):

        #calculate atmospheric densities for given temperature and pressure
        if P is None:
            P = self.P
        if T is None:
            T = self.T # used to be params.planet_temp!
        return  (P)/(KBOLTZ*T)



    def set_mixing_ratios(self):

        # set up mixing ratio array from parameter file inputs

        mixing = self.params.planet_mixing

        X = np.zeros((self.ngas,self.nlayers))
        if self.params.pre_run:
            # if preselector is running
            X += 1e-4
        else:
            #checking if number of mixing ratios = number of gasses
            if len(mixing) is not self.ngas:
                logging.error('Wrong  number of mixing ratios to molecules in parameter file')
                exit()

            # X = np.tile(mixing, [self.nlayers, 1]).transpose()  # @todo check?
            for i in range(self.ngas):
                X[i,:] += float(mixing[i])
        return X


    def readATMfile(self):
    #reads in .atm file
        try:
            out = loadtxt(self.params.in_atm_file)
        except ValueError:
            out = loadtxt(self.params.in_atm_file,comments='*',skiprows=10)
#         OUT[:,2] *= 1000. #converting from km to m

        out = out[argsort(OUT[:,2]),:]

        return out[:,0:3],transpose(out[:,3:])


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
            out = np.zeros((num,self.data.nspecgrid))
            wave = np.transpose(self.data.readfile(path+filelist[0], interpolate=True)[:,0])
        else:
            tmp = self.data.readfile(path+filelist[0], interpolate=False)
            ABSsize = len(tmp[:,0])
            out = np.zeros((num,ABSsize))
            wave = np.transpose(tmp[:,0])

        for i in range(num):
            # print filelist[i]
            out[i,:] = np.transpose(self.data.readfile(path+filelist[i], interpolate=interpolate2grid)[:,1])* 1e-4 #converting cm^2 to m^2

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







