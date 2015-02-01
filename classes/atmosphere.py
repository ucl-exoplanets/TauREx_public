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
AMU   = 1.660538921e-27 #atomic mass to kg

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


        # set mixing ratios
        self.absorbing_gases = self.params.planet_molec
        self.absorbing_gases_X = self.params.planet_mixing
        self.inactive_gases = self.params.planet_inactive_gases
        self.inactive_gases_X = self.params.planet_inactive_gases_X

        # set other atmosphere specific parameters
        self.planet_temp = self.params.planet_temp
        self.planet_mass = self.params.planet_mass

        if self.params.fit_couple_mu:
            self.planet_mu = self.get_coupled_planet_mu()
        else:
            self.planet_mu = self.params.planet_mu

        self.planet_radius = self.params.planet_radius
        self.planet_grav = self.get_surface_gravity()
        self.scaleheight = self.get_scaleheight()
        self.max_pressure = self.params.tp_max_pres

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

        self.num_T_params = 3 #number of free temperature parameters @todo leave it here. used for emission

        if self.params.in_include_rad:
            self.rad        = self.data.rad
        if self.params.in_include_cia:
            self.cia        = self.data.cia

        # get prior bounds. Only valid for emission.
        if self.params.fit_emission:
            self.bounds = self.get_prior_bounds()
            self.fit_params, self.fit_index, self.fit_count = self._setup_parameter_grid(emission=True)



    #class methods

    def get_coupled_planet_mu(self):

        '''
        Get the mean molecular weight (mu) from atmospheric composition
        '''

        mu = 0
        for idx, gasname in enumerate(self.absorbing_gases):
            mu += self.absorbing_gases_X[idx] * self.data.get_molecular_weight(gasname)
        for idx, gasname in enumerate(self.inactive_gases):
            mu += self.inactive_gases_X[idx] * self.data.get_molecular_weight(gasname)

        return mu

    def get_scaleheight(self, T=None, g=None, mu=None):

        if not T:
            T = self.planet_temp
        if not g:
            g = self.planet_grav
        if not mu:
            mu = self.planet_mu

        return (KBOLTZ*T)/(mu*g)

    def get_surface_gravity(self, mass=None, radius=None):

        #calculate surface gravity of planet using Rp and Mp
        # todo surface gravity calculated at 10 mbar radius. Is it ok? Can't really use radius at P0, as it depends on
        # todo the TP profile, which depends on scaleheight, which depends on surface gravity... it's a loop!

        if not mass:
            mass = self.planet_mass
        if not radius:
            radius = self.planet_radius

        return (G * mass) / (radius**2)

    def get_rho(self, T=None, P=None):

        #calculate atmospheric densities for given temperature and pressure
        if P is None:
            P = self.P
        if T is None:
            T = self.T

        return  (P)/(KBOLTZ*T)


    def setup_pta_grid(self):

        #calculate pressure, temperature, altitude grid

        n_scale  = self.params.tp_num_scale # thickness of atmosphere in number of atmospheric scale heights
        max_z = n_scale * self.scaleheight

        #generatinng altitude-pressure array
        pta_arr = np.zeros((self.nlayers,3))
        pta_arr[:,2] = np.linspace(0,max_z,num=self.nlayers) # altitude
        pta_arr[:,0] = self.max_pressure * np.exp(-pta_arr[:,2]/self.scaleheight)
        pta_arr[:,1] = self.planet_temp

        return pta_arr

    def get_gas_fraction(self, gasname):

        # returns the mixing ratio of gasname. The gas can be either an absorber or an inactive gas
        if gasname in self.data.all_absorbing_gases:
            if gasname in self.params.planet_molec:
                index = self.params.planet_molec.index(gasname)
                return self.X[index, 0]

        elif gasname in self.data.all_inactive_gases:
            if gasname in self.params.planet_inactive_gases:
                index = self.params.planet_inactive_gases.index(gasname)
                return self.inactive_gases_X[index]
        return 0

    def get_prior_bounds(self):

        # this is only used in emission now...

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

        # old parameter grid. works for emission though...

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

        mixing = self.absorbing_gases_X

        X = np.zeros((self.ngas,self.nlayers))
        if self.params.pre_run: # @todo FIX POSITION
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





