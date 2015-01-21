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

        #derived values @todo planet gravity depends on radius as well.
        self.scaleheight = self.get_scaleheight(self.params.planet_temp, data.planet_grav, self.params.planet_mu)

        #setting planetary surface gravity
        self.planet_grav = self.get_surface_gravity()

        # build PTA profile and mixing ratio array (pta = pressure, temp, alt; X = mixing ratios of molecules)
        if self.params.in_use_ATMfile: #@ambiguous statement. both tp_var_atm and in_use_ATMfile can be false.

            # reading atmospheric profile from .atm file
            self.nlayers  = self.data.nlayers
            self.ngas     = self.data.ngas
            self.X        = self.data.X
            self.pta      = self.data.pta

        elif self.params.tp_var_atm:

            # atmospheric profile is fitted
            self.nlayers  = int(self.params.tp_atm_levels)
            self.ngas     = int(len(self.params.planet_molec))
            self.X        = self.set_mixing_ratios() # mixing ratios are read from parameter file or set to 1e-4 if preselector = True
            self.pta      = self.setup_pta_grid()

        self.P          = self.pta[:,0] # pressure array
        self.P_bar      = self.P * 1.0e-5 #convert pressure from Pa to bar
        self.T          = self.pta[:,1] # temperature array
        self.z          = self.pta[:,2] # altitude array

        self.rho        = self.get_rho() # assume self.T, self.P

        self.inactive_gases_X = self.params.planet_inactive_gases_X

        if self.params.in_include_rad:
            self.rad        = self.data.rad
        if self.params.in_include_cia:
            self.cia        = self.data.cia

        if self.params.fit_emission or self.params.fit_transmission:
            self.bounds = self.get_prior_bounds()

        if self.params.fit_emission:
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(emission=True)

        if self.params.fit_transmission:
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(transmission=True)

        #self.bulk_composition = self.get_bulk_composition()


    #class methods

    def get_scaleheight(self,T,g,mu):
        return (KBOLTZ*T)/(mu*g)

    def get_surface_gravity(self, mass=None, radius=None):

        #calculate surface gravity of planet using Rp and Mp

        if not mass:
            mass = self.params.planet_mass
        if not radius:
            radius = self.params.planet_radius

        return (G * mass) / (radius**2)

    def get_rho(self, T=None, P=None):

        #calculate atmospheric densities for given temperature and pressure
        if P is None:
            P = self.P
        if T is None:
            T = self.T
        return  (P)/(KBOLTZ*T)



    def setup_pta_grid(self, T=None):

        #calculate pressure, temperature, altitude grid

        max_p    = self.params.tp_max_pres
        n_scale  = self.params.tp_num_scale # thickness of atmosphere in number of atmospheric scale heights
        n_layers = self.nlayers

        # if T is not None:
        if not T:
            T = self.params.planet_temp

        self.planet_grav = self.data.get_surface_gravity(mass=self.params.planet_mass, radius=self.params.planet_radius)
        self.scaleheight = self.get_scaleheight(T,  self.planet_grav, self.params.planet_mu)

        max_z = n_scale * self.scaleheight

        #generatinng altitude-pressure array
        pta_arr = np.zeros((n_layers,3))
        pta_arr[:,2] = np.linspace(0,max_z,num=n_layers) # altitude
        pta_arr[:,0] = max_p * np.exp(-pta_arr[:,2]/self.scaleheight)

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





