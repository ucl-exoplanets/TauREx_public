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
import scipy.special as spe
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
        if self.params.in_use_ATMfile: #@ambiguous statement. both tp_var_atm and in_use_ATMfile can be false.
            # reading atmospheric profile from .atm file
            self.pta      = self.data.pta
            self.X        = self.data.X
            self.nlayers  = self.data.nlayers
            self.ngas     = self.data.ngas
#             self.X        = np.zeros((self.ngas,self.nlayers)) + (1e-5) #setting up initial mixing ratios
            
        elif self.params.tp_var_atm:
            # atmospheric profile is fitted
            self.nlayers  = int(self.params.tp_atm_levels)
            self.ngas     = int(len(self.params.planet_molec))
            self.X        = self.set_mixing_ratios() # mixing ratios are read from parameter file or set to 1e-4 if preselector = True
            self.pta      = self.setup_pta_grid()

    
        

        #setting internal parameters
        self.g          = data.planet_grav #surface gravity
        self.P          = self.pta[:,0] # pressure array
        self.P_bar      = self.P * 1.0e-5 #convert pressure from Pa to bar
        self.T          = self.pta[:,1] # temperature array
        self.z          = self.pta[:,2] # altitude array
        self.rho        = self.get_rho(T=self.T, P=self.P)
        self.sigma_dict = self.data.sigma_dict
        

        #setting optional parameters
        if self.params.in_include_rad:
            self.rad        = self.data.rad
        if self.params.in_include_cia:
            self.cia        = self.data.cia

        if self.params.gen_type == 'emission':
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(emission=True)
            self.TP_type = 'guillot'
            

        if self.params.gen_type == 'transmission':
            self.fit_params, self.fit_index, self.fit_count = self.setup_parameter_grid(transmission=True)
            self.TP_type = 'isothermal'
        
        #setting up TP profile 
        self.TP_profile = self.set_TP_profile()
            
        #setting free parameter prior bounds
        if self.params.fit_emission or self.params.fit_transmission:
                self.bounds = self.get_prior_bounds()

 


    #class methods
    def get_scaleheight(self,T,g,mu):
        '''
        Calculating the atmospheric scale height
        '''
        return (KBOLTZ*T)/(mu*g)


    def setup_pta_grid(self, T=None):
        '''
        Calculate pressure, temperature, altitude grid
        '''
        max_p    = self.params.tp_max_pres
        n_scale  = self.params.tp_num_scale # thickness of atmosphere in number of atmospheric scale heights
        n_layers = self.nlayers

        # get new scale height if T is provided. Otherwise assume default (should be params.planet_temp)

        if T is None:
            T = self.params.planet_temp
        else:
            T = T[0]

        self.planet_grav = self.data.get_surface_gravity(mass=self.params.planet_mass, radius=self.params.planet_radius)
        self.scaleheight = self.get_scaleheight(T,  self.planet_grav, self.params.planet_mu)

        max_z = n_scale * self.scaleheight

        #generatinng altitude-pressure array
        pta_arr = np.zeros((n_layers,3))
        pta_arr[:,2] = np.linspace(0,max_z,num=n_layers) # altitude
        pta_arr[:,0] = max_p * np.exp(-pta_arr[:,2]/self.scaleheight) #pressure
        pta_arr[:,1] = T #temperature
        
        return pta_arr



    def get_prior_bounds(self): 
        '''
        Partially to be moved to parameter file i guess
        '''
        self.X_priors = [self.params.fit_X_low, self.params.fit_X_up] #lower and upper bounds for column density
        self.T_priors = [self.params.planet_temp - self.params.fit_T_low, self.params.planet_temp + self.params.fit_T_up] #lower and upper bounds for temperature


        #setting up bounds for downhill algorithm
        bounds = []
        for i in xrange(self.ngas):
            bounds.append((self.X_priors[0],self.X_priors[1]))
            
        if self.TP_type   == 'isothermal':
            bounds.append((self.T_priors[0],self.T_priors[1])) #isothermal T 
        elif self.TP_type == 'rodgers':
            for i in xrange(self.nlayers):
                bounds.append((self.T_priors[0],self.T_priors[1])) #layer by layer T
        elif self.TP_type == 'guillot':
            bounds.append((self.T_priors[0],self.T_priors[1]))  #T_irr prior
            bounds.append((0.0,1.0))                            #kappa_irr prior
            bounds.append((0.0,1.0))                            #kappa_v1 prior
            bounds.append((0.0,1.0))                            #kappa_v2 prior
            bounds.append((0.0,1.0))                            #alpha prior

        return bounds


    def setup_parameter_grid(self, transmission=False, emission=False):
        '''
        Setting up the parameter grid (variable in length depending on TP-profile model selected)
        
        fit_params    [abundances (ngas),TP-profile (variable)] #TP profile parameters is 1 for isothermal (only T) but varies for differnt models
        fit_count     [no. abundances, no. TP parameters]
        fit_index     index to fitparams for slicing
        '''
        fit_params = []; fit_count = []; fit_index = []

        # setting up mixing ratios for individual gases
        Xmean = np.mean(self.X_priors)
        cgas = 0 #count of gas parameters
        for i in xrange(self.ngas):
            fit_params.append(Xmean)
            cgas += 1
        fit_count.append(cgas)


        #setting up temperature parameters
        ctp = len(self.TP_priors) #count of TP parameters

        if transmission:
            fit_params.append(self.params.planet_temp) #setting temperature to planet equilibrium temperature 

        if emission:
            for i in xrange(ctp):
                fit_params.append(np.mean(self.TP_priors[i])) #setting the mean value from bounds
 
        #generating index list
        cumidx = fit_count[0]
        fit_index.append(cumidx)

        for i in xrange(1,len(fit_count)):
            cumidx += fit_count[i]
            fit_index.append(cumidx)

        return fit_params, fit_index, fit_count


    def get_rho(self, T=None, P=None):
        '''
        Calculate atmospheric densities for given temperature and pressure
        '''
        if P is None:
            P = self.P
        if T is None:
            T = self.T # used to be params.planet_temp!
        return  (P)/(KBOLTZ*T)



    def set_mixing_ratios(self):
        '''
        Set up mixing ratio array from parameter file inputs
        '''
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
   

    #####################
    # Everything below is related to Temperature Pressure Profiles


    def set_TP_profile(self):
        '''
        Decorator supplying TP_profile with correct TP profile function. 
        Only the free parameters will be provided to the function after TP profile 
        is set.   
        '''
        if self.TP_type is 'isothermal':
            profile = self._TP_isothermal
        elif self.TP_type is 'guillot':
            profile = self._TP_guillot2010
        elif self.TP_type is 'rodgers':
            profile = self._TP_rodgers2000
            
        def tpwrap(fit_params):
            return self._TP_profile(profile,fit_params)
        return tpwrap


    def _TP_profile(self, TP_function, fit_params):
        '''
        TP profile decorator. Takes given TP profile function and fitting parameters 
        To generate Temperature, Pressure, Column Density (T, P, X) grids
        
        fit_params depend on TP_profile selected.
        INDEX splits column densities from TP_profile parameters, INDEX = [Chi, TP]
        '''

        fit_index = self.fit_index

        X_params  = fit_params[:fit_index[0]] #splitting to gas parameters
        TP_params = fit_params[fit_index[0]:] #splitting to TP profile parameters

        #Setting up mixing ratio grid. Convert X into Nd arrays
        self.X[:] = 0.0             #using already generated grid
        for i in xrange(self.ngas):
            self.X[i,:] += X_params[i]
            
        #generating TP profile given input TP_function
        T,P = TP_function(TP_params)
        
        return T,P, self.X


        
    def _TP_isothermal(self,TP_params):
        '''
        TP profile for isothermal atmosphere. Follows the old implementation. 
        '''
        self.pta = self.setup_pta_grid(T=TP_params)
        P = self.pta[:,0]; T = self.pta[:,1]
        
        return T,P 
        
        
    def _TP_guillot2010(self,TP_params):
        '''
        TP profile from Guillot 2010, A&A, 520, A27 (equation 49)
        Using modified 2stream approx. from Line et al. 2012, ApJ, 729,93 (equation 19)
        
        Free parameters required: 
            - T_irr    = planet equilibrium temperature (Line fixes this but we keep as free parameter)
            - kappa_ir = mean infra-red opacity
            - kappa_v1 = mean optical opacity one
            - kappa_v2 = mean optical opacity two
            - alpha    = ratio between kappa_v1 and kappa_v2 downwards radiation stream
            
        Fixed parameters:
            - T_int    = internal planetary heat flux (can be largely ignored. line puts it on 200K for hot jupiter)
            - P        = Pressure grid, fixed to self.P
            - g        = surface gravity, fixed to self.g
        '''
        
        #assigning fitting parameters 
        T_irr = TP_params[0], kappa_ir = TP_params[1], kappa_v1 = TP_params[2]; kappa_v2 = TP_params[3]; alpha = TP_params[4]
        
        gamma_1 = kappa_v1/kappa_ir; gamma_2 = kappa_v2/kappa_ir     
        tau = kappa_ir * self.P / self.g
        
        T_int = 200 #@todo internal heat parameter looks suspicious... needs looking at. 
        
        def eta(gamma, tau):
            part1 = 2.0/3.0 + 2.0 / (3.0*gamma) * (1.0 + (gamma*tau/2.0 - 1.0) * np.exp(-1.0 * gamma * tau))
            part2 = 2.0 * gamma / 3.0 * (1.0 - tau**2/2.0) * spe.expn(2,(gamma*tau))
            return part1 + part2
            
        
        T4 = 3.0*T_int**4/4.0 * (2.0/3.0 + tau) + 3.0*T_irr**4/4.0 *(1.0 - alpha) * eta(gamma_1,tau) + 3.0 * T_irr**4/4.0 * alpha * eta(gamma_2,tau)
        T = T4**(1/4)
        return T, self.P
    
    
    def _TP_rodgers2000(self,TP_params):
        '''
        Layer-by-layer temperature - pressure profile retrieval using dampening factor 
        Introduced in Rodgers (2000): Inverse Methods for Atmospheric Sounding (equation 3.26)
        Featured in NEMESIS code (Irwin et al., 2008, J. Quant. Spec., 109, 1136 (equation 19)
        Used in all Barstow et al. papers. 
        
        Free parameters required: 
            - T  = one temperature per layer of Pressure (P)
            
        Fixed parameters: 
            - P  = Pressure grid, fixed to self.P
            - h  = correlation parameter, in scaleheights, Line et al. 2013 sets this to 7, Irwin et al sets this to 1.5
                  may be left as free and Pressure dependent parameter later.
        '''
        
        #assigning parameters 
        T_init = TP_params
        h      = 7.0
        
#         covmatrix = np.zeros((self.nlayers,self.nlayers))
#         for i in range(self.nlayers):
#             for j in range(self.nlayers):
#                 covmatrix[i,j] = (T_init[i] * T_init[j])**(1/2) * np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[j]))/h)  
        
        T = np.zeros((self.nlayers,1))
        for i in xrange(self.nlayers):
            T[i] = np.sum((T_init[i] * T_init[:])**(1/2) * np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[:]))/h))
        
        return T, self.P









##########################
    def get_prior_bounds_old(self):
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

    def setup_parameter_grid_old(self, transmission=False, emission=False):

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



    def TP_profile_old(self, fit_params, T=None, P=None):

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



