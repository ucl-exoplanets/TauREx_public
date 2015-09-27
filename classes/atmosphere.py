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
import sys
import library_general as libgen


# some constants
KBOLTZ = 1.380648813e-23
G = 6.67384e-11
AMU   = 1.660538921e-27 #atomic mass to kg

class atmosphere(base):

    #@profile
    def __init__(self, data, params=None, tp_profile_type=None, covariance=None):

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
        self.inactive_gases = self.params.planet_inactive_gases
        self.inactive_gases_X = self.params.planet_inactive_gases_X
        self.nallgases = len(self.absorbing_gases) + len(self.inactive_gases)

        # set values for parameter files
        self.planet_temp = self.params.planet_temp
        self.max_pressure = self.params.tp_max_pres
        self.planet_radius = self.params.planet_radius
        self.planet_mass = self.params.planet_mass

        self.nlayers  = self.data.nlayers
        self.ngas     = self.data.ngas

        if self.params.in_use_ATMfile:
            # build PTA profile and mixing ratio array from .atm file
            logging.info('Atmospheric PTA grid has been set by .atm file')
            self.X = self.data.X
            self.pta = self.data.pta
            if len(self.X) <> len(self.absorbing_gases):
                logging.error('The number of molecules specified in the .par file is not consistent with'
                             'the mixing ratios in the .atm file %s ' % self.params.in_atm_file)
                sys.exit()
        else:
            self.absorbing_gases_X = self.params.planet_mixing # mixing ratios are read from parameter file
            self.X = self.set_mixing_ratios()

        self.absorbing_gases_X = self.X[:,0] # assume X from lowest layer. This will be used to calculate mu, if couple_mu is True

        # set mu, planet grav and scaleheight
        if self.params.tp_couple_mu:
            self.planet_mu = self.get_coupled_planet_mu()
            logging.info('Coupling mu to composition. Average mu = ' + str(np.average(self.planet_mu)/AMU))

        else:
            self.planet_mu = self.params.planet_mu


        self.planet_grav = self.get_surface_gravity() # planet gravity at the surface
        self.scaleheight = self.get_scaleheight(T = self.planet_temp)

        # determine height of atmosphere using original scale height & planet gravity
        if isinstance(self.scaleheight, float):
            self.max_z = self.params.tp_num_scale * self.scaleheight
        else:
            self.max_z = self.params.tp_num_scale * np.average(self.scaleheight)

        if not self.params.in_use_ATMfile:
            self.pta      = self.setup_pta_grid()

        self.P        = self.pta[:,0] # pressure array
        self.P_bar    = self.P * 1.0e-5 #convert pressure from Pa to bar
        self.T        = self.pta[:,1] # temperature array
        self.z        = self.pta[:,2] # altitude array
        self.rho        = self.get_rho() # assume self.T, self.P

        # reset atmosphere using scale height ang gravity calculated as a function of altitude

        # set cloud parameters
        self.clouds_lower_P = self.params.in_cld_lower_P
        self.clouds_upper_P = self.params.in_cld_upper_P
        self.clouds_m = self.params.in_cld_m
        self.clouds_a = self.params.in_cld_a

        if self.params.in_use_TP_file == True: # todo temporary. Just to feed a given T profile.
            self.T = np.loadtxt(self.params.in_TP_file)[:,0]

         #selecting TP profile to use
        if self.params.gen_type == 'emission':
            if tp_profile_type is None:
                self.TP_type = self.params.tp_type
            else:
                self.TP_type = tp_profile_type

        elif self.params.gen_type == 'transmission':
            if tp_profile_type is None:
                self.TP_type = self.params.tp_type
            else:
                self.TP_type = tp_profile_type

        #setting up TP profile
        self.set_TP_profile()

        #set non-linear sampling grid for hybrid model
        if tp_profile_type == 'hybrid':
            self.hybrid_covmat = covariance
            self.get_TP_sample_grid(covariance,delta=0.05)

    #class methods

    #@profile
    def get_coupled_planet_mu(self, absorbing_gases_X='', inactive_gases_X=''):

        '''
        Get the mean molecular weight (mu) from atmospheric composition
        '''

        # if absorbing_gases_X == '':
        #     absorbing_gases_X = self.absorbing_gases_X

        if inactive_gases_X == '':
            inactive_gases_X = self.inactive_gases_X

        # get mu for each layer
        mu = np.zeros(self.nlayers)

        for i in range(self.nlayers):
            for idx, gasname in enumerate(self.absorbing_gases):
                mu[i] += self.X[idx, i] * self.data.get_molecular_weight(gasname)

            for idx, gasname in enumerate(self.inactive_gases):
                mu[i] += inactive_gases_X[idx] * self.data.get_molecular_weight(gasname)

            logging.debug('Mean molecular weight for layer %i is %.4f' % (i, mu[i]/AMU))

        return mu

    #@profile
    def get_scaleheight(self, T=None, g=None, mu=None):

        if not T:
            T = self.T
        if not g:
            g = self.planet_grav
        if not mu:
            mu = self.planet_mu

        #Tavg = np.average(T)

        return (KBOLTZ*T)/(mu*g)

    #@profile
    def get_surface_gravity(self, mass=None, radius=None):

        #calculate surface gravity of planet using Rp and Mp
        if not mass:
            mass = self.planet_mass
        if not radius:
            radius = self.planet_radius

        return (G * mass) / (radius**2)

    #@profile
    def get_rho(self, T=None, P=None):

        #calculate atmospheric densities for given temperature and pressure
        if P is None:
            P = self.P
        if T is None:
            T = self.T

        return  (P)/(KBOLTZ*T)

    def get_TP_sample_grid(self, covmat, delta=0.02):
        '''
        Calculates non-linear sampling grid for TP profile from provided covariance
        '''

        Pindex = [0]
        for i in range(self.nlayers-1):
            if np.abs(covmat[0,i] - covmat[0,Pindex[-1]]) >= delta:
                Pindex.append(i)
            elif (i - Pindex[-1]) >= 10:
                Pindex.append(i)
        Pindex.append(self.nlayers-1)

        self.P_index = Pindex
        self.P_sample   = self.P[Pindex]

    #@profile
    def setup_pta_grid(self, T=None):

        '''
        Calculate pressure, temperature, altitude grid
        '''

        if T is None:
            T = self.params.planet_temp
        else:
            T = T[0]

        #generatinng altitude-pressure array
        pta_arr = np.zeros((self.nlayers,3))
        pta_arr[:,2] = np.linspace(0, self.max_z, num=self.nlayers) # altitude
        pta_arr[:,0] = self.max_pressure * np.exp(-pta_arr[:,2]/self.scaleheight)
        pta_arr[:,1] = T

        return pta_arr

    #@profile
    def get_gas_fraction(self, gasname):

        # returns the mixing ratio of gasname. The gas can be either an absorber or an inactive gas
        # careful that this returns the gas mixing ratio at the zeroth level.
        if gasname in self.params.all_absorbing_gases:
            if gasname in self.params.planet_molec:
                index = self.params.planet_molec.index(gasname)
                return self.X[index, 0]

        elif gasname in self.params.all_inactive_gases:
            if gasname in self.params.planet_inactive_gases:
                index = self.params.planet_inactive_gases.index(gasname)
                return self.inactive_gases_X[index]
        return 0

    #@profile
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

    def setup_parameter_grid(self):
        '''
        Setting up the parameter grid (variable in length depending on TP-profile model selected)

        fit_params    [abundances (ngas), TP-profile (variable)] #TP profile parameters is 1 for isothermal (only T) but varies for different models
        fit_count     [no. abundances, no. TP parameters]
        fit_index     index to fitparams for slicing
        '''
        fit_params = []; fit_index = []

        # setting up indices for parameter grid
        cpar = len(self.bounds)
        cgas = self.ngas
        ctp  = cpar - cgas
        fit_count = [cgas,ctp]

        #generating index list
        fit_index = [cgas, cpar]

        #setting all initial guesses to the mean of the parameter bounds
        #this can be modified later should it prove to be too awful
        for par in self.bounds:
            fit_params.append(np.mean(par))


        self.fit_params = np.asarray(fit_params)
        self.fit_index  = fit_index
        self.fit_count  = fit_count

#         return np.asarray(fit_params), fit_index, fit_count


    #@profile
    def set_mixing_ratios(self):

        # set up mixing ratio array from parameter file inputs

        mixing = self.absorbing_gases_X

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
        self.X = X

        return X


    def update_atmosphere(self):

        # update surface gravity and scale height
        self.planet_grav = self.get_surface_gravity()
        self.scaleheight = self.get_scaleheight()


        # set altitude array
        if self.params.in_use_ATMfile:
            self.z = self.pta[:,2]
        else:
            self.z = np.linspace(0, self.max_z, num=self.nlayers)

        self.P = self.max_pressure * np.exp(-self.z/self.scaleheight)
        self.P_bar = self.P * 1.0e-5 #convert pressure from Pa to bar
        self.rho = self.get_rho() # update density

    #####################
    # Everything below is related to Temperature Pressure Profiles

    def set_TP_profile(self, profile=None):
        '''
        Decorator supplying TP_profile with correct TP profile function.
        Only the free parameters will be provided to the function after TP profile
        is set.
        '''
        if profile is not None: #implicit or explicit TP_type setting
            self.TP_type = profile

        if self.TP_type == 'isothermal':
            _profile = self._TP_isothermal
        elif self.TP_type == 'guillot':
            _profile = self._TP_guillot2010
        elif self.TP_type == 'rodgers':
            _profile = self._TP_rodgers2000
        elif self.TP_type == 'hybrid':
            _profile = self._TP_hybrid
        elif self.TP_type == '2point':
            _profile = self._TP_2point
        elif self.TP_type == '3point':
            _profile = self._TP_3point
        else:
            logging.error('Invalid TP profile name')

        def tpwrap(fit_params, **kwargs):
            return self._TP_profile(_profile, fit_params, **kwargs)

        self.TP_profile = tpwrap  #setting TP_profile
        self.TP_setup = True #flag that allows TP profile internal things to be run once only. Sets to False after execution

    def _TP_profile(self, TP_function, TP_params, **kwargs):
        '''
        TP profile decorator. Takes given TP profile function and fitting parameters
        To generate Temperature, Pressure, Column Density (T, P, X) grids

        fit_params depend on TP_profile selected.
        INDEX splits column densities from TP_profile parameters, INDEX = [Chi, TP]
        '''

        # X_params  = fit_params[:self.fit_index[0]] #splitting to gas parameters
        # TP_params = TP_params #splitting to TP profile parameters

        #Setting up mixing ratio grid. Convert X into Nd arrays
        # self.X[:] = 0.0             #using already generated grid
        # for i in xrange(self.ngas):
        #     self.X[i,:] += X_params[i]

        #generating TP profile given input TP_function
        T = TP_function(TP_params, **kwargs)

        self.TP_setup = False #running profile setup only on first call

        return T

    def _TP_isothermal(self, TP_params):
        '''
        TP profile for isothermal atmosphere. Follows the old implementation.
        '''

        self.pta = self.setup_pta_grid(T=TP_params)
        T = self.pta[:,1]

        return T

    def _TP_guillot2010(self, TP_params):
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
            - g        = surface gravity, fixed to self.planet_grav
        '''

        #assigning fitting parameters
        T_irr = TP_params[0];
        kappa_ir = TP_params[1];
        kappa_v1 = TP_params[2];
        kappa_v2 = TP_params[3];
        alpha = TP_params[4]

        gamma_1 = kappa_v1/(kappa_ir + 1e-10); gamma_2 = kappa_v2/(kappa_ir + 1e-10)
        tau = kappa_ir * self.P / self.planet_grav

        T_int = 200 #@todo internal heat parameter looks suspicious... needs looking at.

        def eta(gamma, tau):
            part1 = 2.0/3.0 + 2.0 / (3.0*gamma) * (1.0 + (gamma*tau/2.0 - 1.0) * np.exp(-1.0 * gamma * tau))
            part2 = 2.0 * gamma / 3.0 * (1.0 - tau**2/2.0) * spe.expn(2,(gamma*tau))
            return part1 + part2

        T4 = 3.0*T_int**4/4.0 * (2.0/3.0 + tau) + 3.0*T_irr**4/4.0 *(1.0 - alpha) * eta(gamma_1,tau) + 3.0 * T_irr**4/4.0 * alpha * eta(gamma_2,tau)
        T = T4**0.25

        return np.asarray(T)

    def _TP_rodgers2000(self, TP_params, h=None, covmatrix=None):
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
        if h is None:
            h = self.params.tp_corrlength

        #assigning parameters
        T_init = TP_params

#         covmatrix = np.zeros((self.nlayers,self.nlayers))
#         for i in xrange(self.nlayers):
#                 covmatrix[i,:] = np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[:]))/h)

        if covmatrix is None: #if covariance not provided, generate
            if self.TP_setup: #run only once and save
                self.rodgers_covmat = np.zeros((self.nlayers,self.nlayers))
                for i in xrange(self.nlayers):
                    self.rodgers_covmat[i,:] = np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[:]))/h)
        else:
            self.rodgers_covmat = covmatrix


        T = np.zeros((self.nlayers)) #temperature array

        #correlating temperature grid with covariance matrix
        for i in xrange(self.nlayers):
#             covmat  = np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[:]))/h)
            weights = self.rodgers_covmat[i,:] / np.sum(self.rodgers_covmat[i,:])
            T[i]    = np.sum(weights*T_init)

#         pl.figure(3)
#         pl.imshow(self.rodgers_covmat,origin='lower')
#         pl.colorbar()
#         pl.show()

        #sets list of ascii parameter names. This is used in output module to compile parameters.dat
        if self.TP_setup:
            self.TP_params_ascii = []
            for i in xrange(self.nlayers):
                self.TP_params_ascii.append('T_{0}'.format(str(i)))

        return T

    def _TP_hybrid(self, TP_params, h=None):
        '''
        Layer-by-layer temperature pressure profile. It is a hybrid between the _TP_rodgers2000 profile and
        a second (externally calculated) covariance matrix. The external covariance can be calculated
        using a maximum likelihood retrieval of the TP profile before (or similar).
        The hybrid covariance is given by Cov_hyb = (1-alpha) * Cov_TP_rodgers200 + alpha * Cov_external

        Free parameters required:
            - alpha  = weighting parameter between covariance 1 and 2
            - T      = one temperature per layer of Pressure (P)

        Fixed parameters:
            - P  = Pressure grid, fixed to self.P
            - h  = correlation parameter, in scaleheights, Line et al. 2013 sets this to 7, Irwin et al sets this to 1.5
                  may be left as free and Pressure dependent parameter later.
        '''
        if h is None:
            h = self.params.tp_corrlength

        #if self.TP_setup:
        T = np.zeros((self.nlayers)) #temperature array

        #assigning parameters
        alpha  = TP_params[0]

#         alpha = 0.5
#         print 'alpha ',alpha
#         print 'tparams ', TP_params[:10]
        T_init = TP_params[1:]
        covmatrix = self.hybrid_covmat

        #interpolating fitting temperatures to full grid
        T_interp = np.interp(np.log(self.P[::-1]),np.log(self.P_sample[::-1]),T_init[::-1])[::-1]

#         print T_init

        if self.TP_setup: #run only once and save
            self.rodgers_covmat = np.zeros((self.nlayers,self.nlayers))
            for i in xrange(self.nlayers):
                self.rodgers_covmat[i,:] = np.exp(-1.0* np.abs(np.log(self.P[i]/self.P[:]))/h)

        T[:] = 0.0 #temperature array
        cov_hybrid = (1.0-alpha) * self.rodgers_covmat + alpha * covmatrix

#         pl.figure(3)
        #correlating temperature grid with covariance matrix
        for i in xrange(self.nlayers):
            weights = cov_hybrid[i,:] / np.sum(cov_hybrid[i,:])
            T[i]    = np.sum(weights*T_interp)
#             pl.plot(cov_hybrid[i,:])

#         pl.ion()
#         pl.figure(2)
#         pl.clf()
# #         pl.plot(self.T,self.P,'blue')
# #         pl.plot(T_interp,self.P,'k')
#         pl.plot(T_init,self.P_sample,'ro')
#         pl.yscale('log')
# #         xlabel('Temperature')
# #         ylabel('Pressure (Pa)')
#         pl.gca().invert_yaxis()
#         pl.draw()
#
#         pl.ion()
#         pl.figure(3)
#         pl.clf()
#         pl.plot(self.T,self.P,'blue')
# #         pl.plot(T_interp,self.P,'k')
# #         pl.plot(T_init,self.P_sample,'ro')
#         pl.yscale('log')
# #         xlabel('Temperature')
# #         ylabel('Pressure (Pa)')
#         pl.gca().invert_yaxis()
#         pl.draw()



#         pl.figure(3)
#         pl.plot(weights)
#         pl.show()
#         exit()

        return T

    def set_TP_hybrid_covmat(self,covariance):
        '''
        Setting external covariance for _TP_hybrid
        '''
        self.hybrid_covmat = covariance


    def _TP_2point(self,TP_params):
        '''
        Two point TP profile. Simplest possible profile only fitting for the surface temperature and tropopause temp.
        and pressure. Atmosphere above tropopause is isothermal. Temperature is linearly interpolated in log space.
        No inversion possible here.

        Free parameters required:
            - T1 = surface temperature (at 10bar)
            - T2 = temperature at tropopause (given as difference from T1)
            - P1 = pressure at tropopause

        Fixed parameters:
            - Pressure grid (self.P)
        '''

        maxP = np.max(self.P)
        minP = np.min(self.P)

        T_trop = TP_params[0] - TP_params[1]

        P_params = [maxP,TP_params[-1],minP]
        T_params = [TP_params[0],T_trop,T_trop]

        #creating linear T-P profile
        T = np.interp(np.log(self.P[::-1]), np.log(P_params[::-1]), T_params[::-1])

        return T[::-1]


    def _TP_3point(self,TP_params):
        '''
        Same as 2point TP profile but adds one extra point between surface and troposphere

        Free parameters required:
            - T1 = surface temperature (at 10bar)
            - T2 = point 1 temperature (given as difference from T1)
            - P1 = point 1 pressure
            - T3 = point 2 temperature (given as difference from T1)
            - P2 = point 2 pressure

        Fixed parameters
            - Pressure grid (self.P)
        '''

        maxP = np.max(self.P)
        minP = np.min(self.P)

        T_point1 = TP_params[0] - TP_params[1]
        T_point2 = T_point1 - TP_params[2]
        P_params = [maxP,TP_params[-2],TP_params[-1],minP]
        T_params = [TP_params[0],T_point1,T_point2, T_point2]

        #creating linear T-P profile
        T = np.interp(np.log(self.P[::-1]), np.log(P_params[::-1]), T_params[::-1])

        return T[::-1]

