'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Fitting class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''


#loading libraries     
import os
import sys
import shutil
import logging

import numpy as np
from scipy.optimize import minimize

from library_constants import *
from library_general import *

from matplotlib.pylab import *


try:
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

try: 
    import pymc
    pymc_import = True
except:
    pymc_import = False

try:
    from mpi4py import MPI
    MPIimport = True
except ImportError:
    MPIimport = False
    pass

try:
    import library_cythonised_functions as cy_fun
    cythonised = True
except ImportError:
    cythonised = False

cythonised = False # currently disabelling cythonised functions

from library_constants import *
from library_general import *

class fitting(object):

    def __init__(self, forwardmodel):

        # forwardmodel can be either the transmission or the emission object
        logging.info('Initialise object fitting')

        self.params = forwardmodel.params # get params object from forwardmodel
        self.data = forwardmodel.data    # get data object from forwardmodel
        self.atmosphere = forwardmodel.atmosphere    # get atmosphere object from forwardmodel
        self.forwardmodel = forwardmodel

        self.forwardmodel_type = type(forwardmodel).__name__ #used???
        logging.info('Radiative transfer model: %s' % self.forwardmodel_type)

        # MPI support
        if MPIimport:
            self.MPIrank     = MPI.COMM_WORLD.Get_rank()
            self.MPIsize     = MPI.COMM_WORLD.Get_size()
        else:
            self.MPIrank     = 0
            self.MPIsize     = 0

        if not isinstance(self.data.obs_spectrum, (np.ndarray, np.generic)):
            logging.error('You have not provided an input observed spectrum. Cannot continue with fitting.')
            sys.exit()

        # set some folder names, create some folders

        self.dir_mcmc = os.path.join(self.params.out_path, 'MCMC')
        self.dir_multinest = os.path.join(self.params.out_path, 'multinest')

        if self.MPIrank == 0:
            folders = [self.params.out_path, self.dir_mcmc, self.dir_multinest]
            for f in folders:
                if not os.path.isdir(f):
                    logging.info('Create folder %s' % f)
                    os.makedirs(f)

        # set priors, starting values, and fitting parameters
        self.fit_params_names = [] # machine names
        self.fit_params_texlabels = [] # pretty latex labels
        self.fit_params = []
        self.fit_bounds = []

        self.build_fit_params()

        # initialising output tags
        self.DOWN = False
        self.MCMC = False
        self.NEST = False

    def build_fit_params(self):

        # build the fit_params list: input parameters values to all the fitting routines.
        # build the fit_params_names: names of params in fit_params
        # build the fit_params_texlabels: Latex labels of params in fit_params
        # build the fit_bounds list: boundary conditions for fit_params, used to build (uniform) priors
        
 

        if not self.params.gen_ace: # not using the chemically consistent model
            ##########################################################
            # Mixing ratios of active gases either in log/linear space
            if self.params.fit_fit_active_gases:
                for idx, gasname in enumerate(self.params.atm_active_gases):
                    if self.params.fit_mixratio_log: # fit in log space
                        self.fit_params.append(np.log10(self.params.atm_active_gases_mixratios[idx]))
                        self.fit_bounds.append((np.log10(self.params.fit_mixratio_bounds[0]),
                                                np.log10(self.params.fit_mixratio_bounds[1])))
                        self.fit_params_names.append('log_%s' % gasname)
                        self.fit_params_texlabels.append('log(%s)' % tex_gas_label(gasname))

                    else: # fit in linear space
                        self.fit_params.append(self.params.atm_active_gases_mixratios[idx])
                        self.fit_bounds.append((self.params.fit_mixratio_bounds[0],
                                                self.params.fit_mixratio_bounds[1]))
                        self.fit_params_names.append(gasname)
                        self.fit_params_texlabels.append(tex_gas_label(gasname))

            ##########################################################
            # Mixing ratios of N2 either in log/linear space
            if self.params.fit_fit_N2_mixratio:
                    if self.params.fit_mixratio_log: # fit in log space
                        self.fit_params.append(np.log10(self.params.atm_N2_mixratio))
                        self.fit_bounds.append((np.log10(self.params.fit_mixratio_bounds[0]),
                                                np.log10(self.params.fit_mixratio_bounds[1])))
                        self.fit_params_names.append('log_N2')
                        self.fit_params_texlabels.append('log(%s)' % tex_gas_label('N2'))

                    else: # fit in linear space
                        self.fit_params.append(self.params.atm_active_gases_mixratios[idx])
                        self.fit_bounds.append((self.params.fit_mixratio_bounds[0],
                                                self.params.fit_mixratio_bounds[1]))
                        self.fit_params_names.append(gasname)
                        self.fit_params_texlabels.append(tex_gas_label(gasname))

            ##########################################################
            # He/H2 ratio. Fit in log space
            if self.params.fit_fit_He_H2_ratio:
                self.fit_params.append(np.log10(self.params.atm_He_H2_ratio))
                self.fit_bounds.append((np.log10(self.params.fit_He_H2_ratio_bounds[0]),
                                        np.log10(self.params.fit_He_H2_ratio_bounds[1])))
                self.fit_params_names.append('log_H2_He')
                self.fit_params_texlabels.append('log(H$_2$/He)')
                
            ##########################################################
            # Centered-log-ratio transform of gases all gases.
            # NEEDS FIXING OR DEPRECATION!
            #  Consider that we have N2 and H2/He npw... ? Not sure how this fits in.
            # reparametrize the mixing ratios of absorbing and inactive gases using the centered-log-ratio transformation
            # in this case, all mixing ratios of absorbing and inactive gases are always fitted
            # clr = self.get_mixing_ratio_clr()
            # # append all mixing ratios, minus one (as we're using the centered-log-ratio transformation)
            # gasnames = self.forwardmodel.atmosphere.absorbing_gases + self.forwardmodel.atmosphere.inactive_gases
            # for i in range(self.forwardmodel.atmosphere.nallgases - 1):
            #     self.fit_params.append(clr[i])
            #     self.fit_bounds.append((self.params.fit_clr_bounds[0], self.params.fit_clr_bounds[1]))
            #     self.fit_params_names.append('%s_CLR' % gasnames[i])
            #     self.fit_params_texlabels.append('CLR_i' % i)


        #####################################################
        # Chemically consistent model parameters: metallicity and co
        if self.params.gen_ace:
            if self.params.fit_fit_ace_metallicity:
                self.fit_params_names.append('ace_log_metallicity')
                self.fit_params_texlabels.append('log(Metallicity)')
                self.fit_bounds.append((np.log10(self.params.fit_ace_metallicity_bounds[0]),
                                        np.log10(self.params.fit_ace_metallicity_bounds[1])))
                self.fit_params.append(np.log10(self.params.atm_ace_metallicity))
#                 count_ace += 1
            if self.params.fit_fit_ace_co:
                self.fit_params_names.append('ace_co')
                self.fit_params_texlabels.append('C/O')
                self.fit_bounds.append((self.params.fit_ace_co_bounds[0],
                                        self.params.fit_ace_co_bounds[1]))
                self.fit_params.append(self.params.atm_ace_co)
#                 count_ace += 1

        ##########################################################################
        # TP profile parameters
        # Note that each TP_type has to set fit_TP_nparams, the number of fitted parameters in the model
        if self.params.fit_fit_temp:

            T_bounds = (self.params.fit_tp_iso_bounds[0], self.params.fit_tp_iso_bounds[1])
            T_mean = np.mean(((self.params.fit_tp_iso_bounds[0], self.params.fit_tp_iso_bounds[1])))

            if self.forwardmodel.atmosphere.TP_type   == 'isothermal':

                    self.fit_TP_nparams = 1
                    self.fit_params_names.append('T')
                    self.fit_params_texlabels.append('$T$')
                    self.fit_bounds.append((T_bounds[0],T_bounds[1]))
                    self.fit_params.append(T_mean)

            elif self.forwardmodel.atmosphere.TP_type == 'rodgers':

                self.fit_TP_nparams = self.forwardmodel.atmosphere.nlayers

                for i in range(self.forwardmodel.atmosphere.nlayers):
                    self.fit_bounds.append((T_bounds[0],T_bounds[1])) #layer by layer T
                    self.fit_params.append(T_mean)
                    self.fit_params_names.append('T_%i' %i)
                    self.fit_params_texlabels.append('$T_{%i}$' % i)

            elif self.forwardmodel.atmosphere.TP_type == 'hybrid':

                self.fit_TP_nparams = len(self.forwardmodel.atmosphere.P_index) + 1

                self.fit_params_names.append('alpha')
                self.fit_params_texlabels.append('$\\alpha$')
                self.fit_bounds.append((self.params.fit_hybrid_alpha_l,self.params.fit_hybrid_alpha_h)) #alpha parameter
                self.fit_params.append(np.mean((self.params.fit_hybrid_alpha_l,self.params.fit_hybrid_alpha_h)))

                for i in range(len(self.forwardmodel.atmosphere.P_index)):
                    self.fit_params_names.append('T_%i' % i)
                    self.fit_params_texlabels.append('$T_{%i}$' % i)
                    self.fit_bounds.append((T_bounds[0],T_bounds[1])) #layer by layer T
                    self.fit_params.append(T_mean)

            elif self.forwardmodel.atmosphere.TP_type == 'guillot':

                self.fit_TP_nparams = 5

                self.fit_params_names.append('T_irr')
                self.fit_params_texlabels.append('$T_\mathrm{irr}$')
                self.fit_bounds.append((self.params.fit_tp_guillot_T_irr_bounds[0], self.params.fit_tp_guillot_T_irr_bounds[1]))
                self.fit_params.append(np.mean((self.params.fit_tp_guillot_T_irr_bounds[0], self.params.fit_tp_guillot_T_irr_bounds[1])))

                self.fit_params_names.append('kappa_ir')
                self.fit_params_texlabels.append('$k_\mathrm{ir}$')
                self.fit_bounds.append((np.log10(self.params.fit_tp_guillot_kappa_ir_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_ir_bounds[1])))
                self.fit_params.append(np.mean((np.log10(self.params.fit_tp_guillot_kappa_ir_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_ir_bounds[1]))))

                self.fit_params_names.append('kappa_v1') #
                self.fit_params_texlabels.append('$k_\mathrm{1}$')
                self.fit_bounds.append((np.log10(self.params.fit_tp_guillot_kappa_v1_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_v1_bounds[1])))
                self.fit_params.append(np.mean((np.log10(self.params.fit_tp_guillot_kappa_v1_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_v1_bounds[1]))))

                self.fit_params_names.append('kappa_v2')
                self.fit_params_texlabels.append('$k_\mathrm{2}$')
                self.fit_bounds.append((np.log10(self.params.fit_tp_guillot_kappa_v2_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_v2_bounds[1])))
                self.fit_params.append(np.mean((np.log10(self.params.fit_tp_guillot_kappa_v2_bounds[0]), np.log10(self.params.fit_tp_guillot_kappa_v2_bounds[1]))))

                self.fit_params_names.append('alpha')
                self.fit_params_texlabels.append('$\\alpha$')
                self.fit_bounds.append((self.params.fit_tp_guillot_alpha_bounds[0], self.params.fit_tp_guillot_alpha_bounds[1]))
                self.fit_params.append(np.mean((self.params.fit_tp_guillot_alpha_bounds[0], self.params.fit_tp_guillot_alpha_bounds[1])))

            elif self.forwardmodel.atmosphere.TP_type == '2point':

                self.fit_TP_nparams = 3

                self.fit_params_names.append('T_surf') #surface layer T
                self.fit_params_texlabels.append('$T_mathrm{surf}$')
                self.fit_bounds.append((T_bounds[0],T_bounds[1]))
                self.fit_params.append(T_mean)

                self.fit_params_names.append('T_trop') #troposphere layer T difference (T_surface- Tdiff) = T_trop
                self.fit_params_texlabels.append('$T_mathrm{trop}$')
                self.fit_bounds.append((0.0,1000.0))
                self.fit_params.append(np.mean((0.0,1000.0)))

                self.fit_params_names.append('P_trop') #troposphere pressure (Pa) #@todo careful with this needs to move somewhere else
                self.fit_params_texlabels.append('$P_mathrm{trop}$')
                self.fit_bounds.append((1.0,1e5))
                self.fit_params.append(np.mean((1.0,1e5)))

            elif self.forwardmodel.atmosphere.TP_type == '3point':

                self.fit_TP_nparams = 5

                self.fit_params_names.append('T_surf')
                self.fit_params_texlabels.append('$T_mathrm{surf}$')
                self.fit_bounds.append((T_bounds[0],T_bounds[1])) #surface layer T
                self.fit_params.append(T_mean)

                self.fit_params_names.append('T_point1') #point1 T difference (T_surface- Tdiff) = T_point1
                self.fit_params_texlabels.append('$T_1$')
                self.fit_bounds.append((0.0,500.0))
                self.fit_params.append(np.mean((0.0,500.0)))

                self.fit_params_names.append('T_point2') #point2 T difference (T_point1- Tdiff) = T_point2
                self.fit_params_texlabels.append('$T_2$')
                self.fit_bounds.append((0.0,500.0))
                self.fit_params.append(np.mean((0.0,500.0)))

                self.fit_params_names.append('P_point1')  #point1 pressure (Pa) #@todo careful with this needs to move somewhere else
                self.fit_params_texlabels.append('$P_1$')
                self.fit_bounds.append((1.0,1e5))
                self.fit_params.append(np.mean((1.0,1e5)))

                self.fit_params_names.append('P_point2') #point2 pressure (Pa) #@todo careful with this needs to move somewhere else
                self.fit_params_texlabels.append('$P_2$')
                self.fit_bounds.append((1.0,1e5))
                self.fit_params.append(np.mean((1.0,1e5)))
                

        else: # not fitting for the TP profile
            self.fit_TP_nparams = 0

        # ##########################################################################
        # DEPRECATED
        # # mean molecular weight. Only if we are not coupling mu to the mixing ratios
        # if not self.params.fit_couple_mu:
        #     if self.params.fit_fit_mu:
        #         self.fit_params_names.append('mu')
        #         self.fit_params_texlabels.append('$\mu$')
        #         self.fit_params.append(self.forwardmodel.atmosphere.mu_profile[0]/AMU) # in AMU
        #         self.fit_bounds.append((self.params.fit_mu_bounds[0], self.params.fit_mu_bounds[1])) # in AMU

        ##########################################################################
        # Fit for 10 bar radius
        if self.params.fit_fit_radius:
            self.fit_params_names.append('radius')
            self.fit_params_texlabels.append('$R_p$')
            self.fit_params.append(self.forwardmodel.atmosphere.planet_radius/RJUP)
            if self.params.fit_radius_bounds_factor:
                # If radius_bounds_factor != None/False, radius bounds are defined as +/- a factor of the input 10 bar radius.
                radius_lower = self.forwardmodel.atmosphere.planet_radius/RJUP -\
                               self.forwardmodel.atmosphere.planet_radius/RJUP * self.params.fit_radius_bounds_factor
                radius_upper = self.forwardmodel.atmosphere.planet_radius/RJUP +\
                               self.forwardmodel.atmosphere.planet_radius/RJUP * self.params.fit_radius_bounds_factor
                self.fit_bounds.append((radius_lower,  radius_upper)) # in RJUP

            else:
                # radius bounds are defined explicitly (radius_bounds)
                self.fit_bounds.append((self.params.fit_radius_bounds[0],  self.params.fit_radius_bounds[1])) # in RJUP

        # ##########################################################################
        # # surface pressure
        # DEPRECATED. USE CLOUD PRESSURE INSTEAD!
        # if self.params.fit_fit_P0:
        #     self.fit_params_names.append('P0')
        #     self.fit_params_texlabels.append('$P_\mathrm{max}$')
        #     self.fit_params.append(self.params.atm_max_pres)
        #     self.fit_bounds.append((np.log10(self.params.fit_P0_bounds[0]), np.log10(self.params.fit_P0_bounds[1]))) # in log[Pascal]

        ##########################################################################
        # Cloud pressure
        if self.params.fit_fit_clouds_pressure:
            self.fit_params_names.append('clouds_pressure')
            self.fit_params_texlabels.append('$P_\mathrm{clouds}$')
            self.fit_params.append(np.mean((np.log10(self.params.fit_clouds_pressure_bounds[0]),
                                            np.log10(self.params.fit_clouds_pressure_bounds[1]))))
            self.fit_bounds.append((np.log10(self.params.fit_clouds_pressure_bounds[0]),
                                    np.log10(self.params.fit_clouds_pressure_bounds[1])))

        # define total number of parameters to be fitted
        self.fit_nparams = len(self.fit_params)
        
        #Mie scattering. Replaces clouds and Rayleigh
        if self.params.fit_fit_mie:
            
            self.fit_params_names.append('log_cloud_mixing')
            self.fit_params_texlabels.append('log($\chi_\mathrm{clouds}$)')
            self.fit_params.append(np.mean((np.log10(self.params.fit_mie_f_bounds[0]),
                                                np.log10(self.params.fit_mie_f_bounds[1]))))
            self.fit_bounds.append((np.log10(self.params.fit_mie_f_bounds[0]),
                                        np.log10(self.params.fit_mie_f_bounds[1])))
            
            
            
            if self.params.fit_fit_mie_radius:
                self.fit_params_names.append('log_clouds_particle_size')
                self.fit_params_texlabels.append('log($R_\mathrm{clouds}$)')
                self.fit_params.append(np.mean((np.log10(self.params.fit_mie_r_bounds[0]),
                                                np.log10(self.params.fit_mie_r_bounds[1]))))
                self.fit_bounds.append((np.log10(self.params.fit_mie_r_bounds[0]),
                                        np.log10(self.params.fit_mie_r_bounds[1])))
            
            if self.params.fit_fit_mie_composition and self.params.atm_mie_type == 'lee':
                self.fit_params_names.append('clouds_composition')
                self.fit_params_texlabels.append('$Q_\mathrm{clouds}$')
                self.fit_params.append(np.mean((self.params.fit_mie_q_bounds[0],
                                                    self.params.fit_mie_q_bounds[1])))
                self.fit_bounds.append((self.params.fit_mie_q_bounds[0],
                                            self.params.fit_mie_q_bounds[1]))
                
            if self.params.fit_fit_mie_cloud_topP:
                self.fit_params_names.append('log_clouds_topP')
                self.fit_params_texlabels.append('$log(P_\mathrm{top})$')
                self.fit_params.append(np.mean((np.log10(self.params.fit_mie_topP_bounds[0]),
                                                    np.log10(self.params.fit_mie_topP_bounds[1]))))
                self.fit_bounds.append((np.log10(self.params.fit_mie_topP_bounds[0]),
                                            np.log10(self.params.fit_mie_topP_bounds[1])))
                
            if self.params.fit_fit_mie_cloud_bottomP:
                self.fit_params_names.append('log_clouds_bottomP')
                self.fit_params_texlabels.append('$log(P_\mathrm{bottom})$')
                self.fit_params.append(np.mean((np.log10(self.params.fit_mie_bottomP_bounds[0]),
                                                    np.log10(self.params.fit_mie_bottomP_bounds[1]))))
                self.fit_bounds.append((np.log10(self.params.fit_mie_bottomP_bounds[0]),
                                            np.log10(self.params.fit_mie_bottomP_bounds[1])))
                                                
        logging.info('Dimensionality: %i' % len(self.fit_params_names))
        logging.info('Fitted parameters name: %s' % self.fit_params_names)
        logging.info('Fitted parameters value: %s' % self.fit_params)
        logging.info('Fitted parameters bound: %s' % self.fit_bounds)
        
        
    #@profile
    def update_atmospheric_parameters(self, fit_params):

        ''' update atmospheric parameters with new values from fit_params
        fit_params contains list of parameters as defined in build_fit_params
        the order in which parameters are read from fit_params must be the same as how these were
        defined in build_fit_params '''

        if type(fit_params).__name__ == 'ndarray':
            fit_params = fit_params.tolist()

        ###########################################################################
        # Parameters in fit_params are stored in this order:
        #
        # Param name                     Number of parameters              Condition
        # -------------------------------------------------------------------------------------------------------------
        # log_GASNAME or GASNAME         atmosphere.nactivegases           params.fit_fit_active_gases = True
        # log_N2 or N2                   1                                 params.fit_fit_N2_mixratio = True
        # log_H2_He                      1                                 params.fit_fit_He_H2_ratio = True
        # ace_log_metallicity            1                                 params.gen_ace = True
        #                                                                  (disable fit of previous params)
        # ace_co                         1                                 ...
        # temperature parameters         Depending on params.atm_tp_type   params.fit_fit_temp = True
        # radius                         1                                 params.fit_fit_radius = True
        # clouds_pressure                1                                 params.fit_fit_clouds_pressure = True
        # mie_scattering                 3 - 4                             params.fit_fit_mie = True


        count = 0 # used to iterate over fit_params[count]

        if not self.params.gen_ace: # not using the chemically consistent model

            ##########################################################
            # Mixing ratios of active gases either in log/linear space
            if self.params.fit_fit_active_gases:
                for idx, gasname in enumerate(self.params.atm_active_gases):
                    if self.params.fit_mixratio_log: # fit in log space
                        self.forwardmodel.atmosphere.active_mixratio_profile[idx, :] = np.power(10, fit_params[count])
                    else:
                        self.forwardmodel.atmosphere.active_mixratio_profile[idx, :] = fit_params[count]
                    count += 1

            ##########################################################
            # Mixing ratios of N2 either in log/linear space
            if self.params.fit_fit_N2_mixratio:
                # note that N2 is idx = 2 in inactive_mixratio_profile
                if self.params.fit_mixratio_log: # fit in log space
                    self.forwardmodel.atmosphere.inactive_mixratio_profile[2, :] = np.power(10, fit_params[count])
                else:
                    self.forwardmodel.atmosphere.inactive_mixratio_profile[2, :] = fit_params[count]
                count += 1

            ##########################################################
            # H2/He ratio. Fit in log space
            if self.params.fit_fit_He_H2_ratio:

                He_H2_ratio = np.power(10, fit_params[count])
                #He_H2_ratio = fit_params[count]

                # first get the sum of the mixing ratio of all active gases
                if len(self.atmosphere.active_gases) > 1:
                    active_mixratio_sum = np.sum(self.atmosphere.active_mixratio_profile, axis = 0)
                else:
                    active_mixratio_sum = np.copy(self.atmosphere.active_mixratio_profile)

                # add the N2 mixing ratio profile to this sum
                active_mixratio_sum += self.atmosphere.inactive_mixratio_profile[2, :]

                # error if this sum is larger than 1 (i.e. 100%)
                # if np.any(active_mixratio_sum > 1.):
                #     # find a way to strongly penalize this in the fitting?

                # The remainder of the atmosphere is made of a mixture of He/H2, with ratio given by the He_H2_ratio
                mixratio_remainder = 1. - active_mixratio_sum
                self.atmosphere.inactive_mixratio_profile[0, :] = mixratio_remainder/(1. + He_H2_ratio) # H2
                self.atmosphere.inactive_mixratio_profile[1, :] =  He_H2_ratio * \
                                                                   self.atmosphere.inactive_mixratio_profile[0, :]

                count += 1

        ####################################################################################
        # Mixing ratios of all gases are transformed using a centered-log-ratio transformation.
        # * NOT WORKING
        # if not self.params.gen_ace: # only if we are not running the chemically consistent model
        #     if self.params.fit_clr_trans: # running centered-log-ratio transformation applied
        #         # convert centered-log-ratio mixing ratios back to simplex
        #         # build list of log-ratios
        #         clr = fit_params[:self.forwardmodel.atmosphere.nallgases-1]
        #         clr.append(-sum(clr)) # append last log-ratio. This is not fitted, but derived as sum(log-ratios) = 0
        #         # convert log-ratios to simplex
        #         clr_inv = self.get_mixing_ratio_inv_clr(clr) # @todo could be improved
        #         # set mixing ratios of absorbing and inactive gases
        #         j = 0
        #         for idx, gasname in enumerate(self.params.atm_active_gases):
        #             self.forwardmodel.atmosphere.atm_active_gases_mixratios[idx,:] = clr_inv[j]
        #             j += 1
        #         for idx, gasname in enumerate(self.params.atm_inactive_gases):
        #             self.forwardmodel.atmosphere.inactive_gases_X[idx] = clr_inv[j]
        #             j += 1
        #         count = self.forwardmodel.atmosphere.nallgases - 1


        ####################################################################################
        # Chemically consistent model parameters. Gas abundances are set by the chemically consistent model ACE.
        if self.params.gen_ace:
            if self.params.fit_fit_ace_metallicity: # metallicity
                self.atmosphere.ace_metallicity = np.power(10, fit_params[count])
                count += 1
            if self.params.fit_fit_ace_co: # c/o ratio
                self.atmosphere.ace_co = fit_params[count]
                count += 1

        ####################################################################################
        # Temperature-pressure profile
        # get TP profile fitted parameters. Number of parameter is profile dependent, and defined by self.fit_TP_nparams
        if self.fit_TP_nparams > 0:
            TP_params = fit_params[count:count+self.fit_TP_nparams]
            self.forwardmodel.atmosphere.temperature_profile = self.forwardmodel.atmosphere.TP_profile(fit_params=TP_params)
            count += self.fit_TP_nparams

        # #####################################################
        # # Mean molecular weight.
        # # DEPRECATED. MMW IS SET BEFORE
        # #    If coupling, then we just derive mu from the mixing ratios.
        # #    If we're fitting, get it from fit_params
        # if self.params.fit_couple_mu:
        #     # mu is not fitted, but coupled to the atmospheric composition
        #     self.forwardmodel.atmosphere.mu_profile = self.forwardmodel.atmosphere.get_coupled_mu_profile()
        # else:
        #     if self.params.fit_fit_mu:
        #         # mu is fitted for the remainder of the atmosphere (i.e. the mu of the atmosphere excluding the active gases)
        #         # the fitted mu is not the mmw of the *entire* atmosphere, but of the remainder of the atmosphere (1 - sum(active gases mixratio) )
        #         mw_active_gases = 0
        #         mixratio_active_gases = 0
        #         # calculate contribution to total mu from active gases
        #         for idx, gasname in enumerate(self.params.atm_active_gases):
        #             mixratio_active_gases += self.forwardmodel.atmosphere.active_mixratio_profile[idx, 0]
        #             mw_active_gases += self.forwardmodel.atmosphere.active_mixratio_profile[idx, 0] * get_molecular_weight(gasname) # assume mu from first layer
        #         # set global mu of the entire atmosphere
        #         self.forwardmodel.atmosphere.mu_profile[:]= mw_active_gases + (1.-mixratio_active_gases)*fit_params[count]*AMU
        #         # see the parameter Fitting>inactive_mu_rescale in default.par for explanation
        #         if self.params.fit_inactive_mu_rescale:
        #             mw_0 = get_molecular_weight(self.params.atm_inactive_gases[0])
        #             mw_1 = get_molecular_weight(self.params.atm_inactive_gases[1])
        #             mixratio_1 = (fit_params[count]*AMU - mw_1)/(mw_0 - mw_1)
        #             mixratio_2 = 1. - mixratio_1
        #             mixratio_remainder = 1. - np.sum(self.forwardmodel.atmosphere.active_mixratio_profile[:,0], axis=0)
        #             inactive_mixratio_profile = np.zeros((self.forwardmodel.atmosphere.ninactivegases, self.atmosphere.nlayers))
        #             self.forwardmodel.atmosphere.inactive_mixratio_profile[0, :] = mixratio_remainder*mixratio_1
        #             self.forwardmodel.atmosphere.inactive_mixratio_profile[1, :] = mixratio_remainder*mixratio_2
        #
        #         count += 1

        ####################################################################################
        # Radius
        if self.params.fit_fit_radius:
            self.forwardmodel.atmosphere.planet_radius = fit_params[count]*RJUP
            count += 1

        # #####################################################
        # # Surface pressure
        # DEPRECATED. Use clouds_pressure
        # if self.params.fit_fit_P0:
        #     self.forwardmodel.atmosphere.max_pressure = np.power(10, fit_params[count])
        #     count += 1

        ####################################################################################
        # Clouds pressure
        if self.params.fit_fit_clouds_pressure:
            self.forwardmodel.atmosphere.clouds_pressure = np.power(10, fit_params[count])
            count += 1
            
        ####################################################################################
        # Mie scattering
        if self.params.fit_fit_mie:
            self.forwardmodel.atmosphere.mie_f = np.power(10,fit_params[count])
            count += 1
            if self.params.fit_fit_mie_radius:
                self.forwardmodel.atmosphere.mie_r = np.power(10,fit_params[count])
                count += 1
            
            if self.params.fit_fit_mie_composition:
                self.forwardmodel.atmosphere.mie_q = fit_params[count]
                count += 1
            
            if self.params.fit_fit_mie_cloud_topP:
                self.forwardmodel.atmosphere.mie_topP = np.power(10,fit_params[count])
                count += 1
            if self.params.fit_fit_mie_cloud_bottomP:
                self.forwardmodel.atmosphere.mie_topP = np.power(10,fit_params[count])
                count += 1

        # if self.params.fit_fit_P0:
        # DEPRECTED
        #     self.forwardmodel.atmosphere.pressure_profile = self.forwardmodel.atmosphere.get_pressure_profile()

        ####################################################################################
        # END. All parameters have been extracted from fit_params
        ####################################################################################

        # Update the state of the atmosphere

        self.atmosphere.set_mu_profile() # update the planet mmw

        if self.fit_TP_nparams > 0:  # density profile changes only if the temperature changes
            self.forwardmodel.atmosphere.set_density_profile()

        if not self.params.gen_ace:
            # update altitude, gravity and scale height profile
            self.forwardmodel.atmosphere.set_altitude_gravity_scaleheight_profile()
        else:
            # if running the chiemical consistent model, set the corresponding parameters (this will use
            # atmosphere.ace_metallicity and atmosphere.ace_co set above
            self.forwardmodel.atmosphere.set_ace_params()
            # Note that if gen_ace is True, set_altitude_gravity_scaleheight_profile is called from the
            # transmission/emission object
        if self.params.fit_fit_mie: #updates mie scattering slope
            self.forwardmodel.atmosphere.get_mie_opacities()


    #@profile
    def chisq_trans(self, fit_params, data, datastd):

        # update atmospheric parameters with fit_params values
        self.update_atmospheric_parameters(fit_params)
        
        # get forward model and bin
        model_out = self.forwardmodel.model()

        if self.params.in_opacity_method in ['xsec_sampled', 'xsec_lowres', 'xsec']:
            # bin if using sampled cross sections
            model = [model_out[self.atmosphere.int_bingrididx == i].mean() for i in range(1, self.atmosphere.int_nbingrid+1)]
            #model_binned = cy_fun.runtime_bin_spectrum(model, self.atmosphere.int_bingrididx, self.atmosphere.int_nbingrid)
        elif self.params.in_opacity_method in ['ktab', 'ktable', 'ktables']:
            # interpolate if using ktables
            model = np.interp(self.data.obs_wngrid, self.atmosphere.int_wngrid, model_out)

        # get chi2
        res = ((data - model) / datastd)
        
        res = np.nansum(res*res)
        if res == 0:
            res = np.nan

           
#         ion()
#         figure(1)
#         clf()
#         plot(self.forwardmodel.atmosphere.temperature_profile, self.forwardmodel.atmosphere.pressure_profile)
#         gca().invert_yaxis()
#         xlabel('Temperature')
#         ylabel('Pressure (Pa)')
#         yscale('log')
#         draw()
#         figure(2)
#         clf()
#    
#         ion()
#         clf()
#         errorbar(self.data.obs_spectrum[:,0],self.data.obs_spectrum[:,1],self.data.obs_spectrum[:,2])
#         plot(self.data.obs_spectrum[:,0], model)
#         xlabel('Wavelength (micron)')
#         ylabel('Transit depth')
#         xscale('log')
#         xlim((min(self.data.obs_spectrum[:,0]), max(self.data.obs_spectrum[:,0])))
#         draw()
#         pause(0.0001)
# # # # # 
#         print('res=%.1f - T=%.1f, mu=%.2f, R=%.4f,' % (res, self.forwardmodel.atmosphere.temperature_profile[0], \
#             self.forwardmodel.atmosphere.mu_profile[0]/AMU, \
#             self.forwardmodel.atmosphere.planet_radius/RJUP), \
#             fit_params) #fit_params

        return res

    ###############################################################################
    #simplex downhill algorithm

    def downhill_fit(self):
    # fitting using downhill algorithm. 
    # @todo makes this proper MLE by replacing chi-squared with likelihood function

        logging.info('Fit data using %s minimisation' % self.params.downhill_type)

        data = self.data.obs_spectrum[:,1] # observed data
        datastd = self.data.obs_spectrum[:,2] # data error

        fit_output = minimize(fun=self.chisq_trans,
                              x0=self.fit_params,
                              args=(data,datastd),
                              method=self.params.downhill_type,
                              bounds=(self.fit_bounds))

        logging.info('Saving the downhill minimization results')

        self.DOWN_fit_total = fit_output
        self.DOWN_fit_output = fit_output['x']
        self.DOWN = True


    ###############################################################################
    #Markov Chain Monte Carlo algorithm. OUTPUT not currently supported

    def mcmc_fit(self):
    # fitting using adaptive Markov Chain Monte Carlo

        logging.info('Start MCMC fit')

        if self.DOWN:
            self.fit_params_init = self.DOWN_fit_output
        else:
            self.fit_params_init = self.fit_params

        data = self.data.obs_spectrum[:,1] #observed data
        datastd = self.data.obs_spectrum[:,2] #data error

        # setting prior distributions
        # master thread (MPIrank =0) will start from ideal solution (from downhill fitting)
        # slave threads (MPIrank != 0) will start from randomised points.

        priors = empty(size(self.fit_params_init), dtype=object)

        if self.MPIrank == 0:
            # setting up main thread. Use downhill FIT as starting points
            for i in range(self.fit_nparams):
                    priors[i] = pymc.Uniform('PFIT_%i' % (i), self.fit_bounds[i][0], self.fit_bounds[i][1],
                                             value=self.fit_params_init[i])  # uniform prior
        else:
            #setting up other threads (if exist). Their initial starting positions will be randomly perturbed
            for i in range(self.fit_nparams):
                param_range = (self.fit_bounds[i][1] - self.fit_bounds[i][0]) / 5.0 #range of parameter over which to perturb starting position
                param_mean  = np.mean(self.fit_bounds[i])
                param_rand  = random.uniform(low=param_mean-param_range,high=param_mean+param_range) #random parameter start
                priors[i] = pymc.Uniform('PFIT_%i' % (i), self.fit_bounds[i][0],self.fit_bounds[i][1],value=param_rand)  # uniform prior

        #setting up data error prior if specified
        if self.params.mcmc_update_std:
            #uniform prior on data standard deviation
            std_dev = pymc.Uniform('datastd', 0.0, 2.0*max(datastd), value=datastd, size=len(datastd))
        else:
            std_dev = pymc.Uniform('datastd', 0.0, 2.0*max(datastd), value=datastd, observed=True, size=len(datastd))

        def mcmc_loglikelihood(value, fit_params, datastd, data):
            # log-likelihood function. Needs to be initialised directly since CYTHON does not like PYMC decorators
            # @todo need to cq ast from numpy object array to float array. Slow but hstack is  slower.
            fit_params_container = zeros((self.fit_nparams))
            for i in range(self.fit_nparams):
                fit_params_container[i] = fit_params[i]
            # @todo params_container should be equal to PFIT? I think so... Maybe not if we fix some values
            chi_t = self.chisq_trans(fit_params_container, data, datastd) #calculate chisq
            llterms =   (-len(data)/2.0)*log(pi) - log(mean(datastd)) - 0.5* chi_t
            return llterms

        mcmc_logp = pymc.Stochastic(logp = mcmc_loglikelihood,
                                    doc = 'The switchpoint for mcmc loglikelihood.',
                                    name = 'switchpoint',
                                    parents = {'fit_params': priors, 'datastd': datastd, 'data': data},
                                    #random = switchpoint_rand,
                                    trace = True,
                                    value = self.fit_params_init,
                                    dtype=int,
                                    rseed = 1.,
                                    observed = True,
                                    cache_depth = 2,
                                    plot=False,
                                    verbose = 0)

        # set output folder
        dir_mcmc_thread = os.path.join(self.dir_mcmc, 'thread_%i' % self.MPIrank)

        # remove old files
        if os.path.isdir(dir_mcmc_thread):
            logging.debug('Remove folder %s' % dir_mcmc_thread)
            shutil.rmtree(dir_mcmc_thread)

        # executing MCMC sampling
        if self.params.mcmc_verbose:
            verbose = 1
        else:
            verbose = 0

        # build the model
        R = pymc.MCMC((priors, mcmc_logp),
                      verbose=verbose,
                      db='txt',
                      dbname=dir_mcmc_thread)

        # populate and run it
        R.sample(iter=self.params.mcmc_iter,
                 burn=self.params.mcmc_burn,
                 thin=self.params.mcmc_thin,
                 progress_bar=self.params.mcmc_progressbar,
                 verbose=verbose)

        # todo Save similar output: NEST_out
        self.MCMC_R = R
        self.MCMC = True

    ###############################################################################
    #Nested Sampling algorithm

    def multinest_fit(self, resume=None):

        logging.info('Start MULTINEST fit')

        if resume is None:
            resume = self.params.nest_resume

        if resume:
            # Resuming from previous run. Restart from scratch if number of live points has changed.
            filename = os.path.join(self.dir_multinest, '1-live.points')
            if os.path.isfile(filename):
                livepoints = sum(1 for line in open(filename))
                if livepoints != self.params.nest_nlive:
                    # todo this should only be done if MPIrank = 0
                    logging.warning('Cannot resume previous MULTINEST run, the number of live points has changed')
                    logging.warning('Removing previous MULTINEST chains')
                    for file in os.listdir(self.dir_multinest):
                        file_path = os.path.join(self.dir_multinest, file)
                        try:
                            if os.path.isfile(file_path):
                                os.unlink(file_path)
                        except:
                            logging.error('Cannot remove files in %s' % self.dir_multinest)

        def multinest_loglike(cube, ndim, nparams):
            # log-likelihood function called by multinest
            fit_params_container = asarray([cube[i] for i in range(len(self.fit_params))])
            chi_t = self.chisq_trans(fit_params_container, data, datastd)
            loglike = (-1.)*np.sum(np.log(datastd*np.sqrt(2*np.pi))) - 0.5 * chi_t
            return loglike

        def multinest_uniform_prior(cube, ndim, nparams):
            # prior distributions called by multinest. Implements a uniform prior
            # converting parameters from normalised grid to uniform prior
            for i in range(len(self.fit_params)):
                cube[i] = (cube[i] * (self.fit_bounds[i][1]-self.fit_bounds[i][0])) + self.fit_bounds[i][0]

        data = self.data.obs_spectrum[:,1] # observed data
        datastd = self.data.obs_spectrum[:,2] # data error
        datastd_mean = mean(datastd)
        ndim = len(self.fit_params)

        logging.info('Multinest output dir is %s' % self.dir_multinest)

        pymultinest.run(LogLikelihood=multinest_loglike,
                        Prior=multinest_uniform_prior,
                        n_dims=len(self.fit_params),
                        multimodal=self.params.nest_multimodes,
                        max_modes=self.params.nest_max_modes,
                        outputfiles_basename=os.path.join(self.dir_multinest, '1-'),
                        const_efficiency_mode = self.params.nest_const_eff,
                        importance_nested_sampling = self.params.nest_imp_sampling,
                        resume = resume,
                        verbose = self.params.nest_verbose,
                        sampling_efficiency = self.params.nest_samp_eff,
                        evidence_tolerance = self.params.nest_ev_tol,
                        mode_tolerance = self.params.nest_mode_tol,
                        n_live_points = self.params.nest_nlive,
                        max_iter= self.params.nest_max_iter,
                        init_MPI=False)

        # wait for all threads to synchronise
        if MPIimport:
            MPI.COMM_WORLD.Barrier()

        self.NEST = True

    # def get_mixing_ratio_clr(self):
    #
    #     # get the geometrical mean
    #     sumlog = 0
    #     n = len(self.forwardmodel.atmosphere.absorbing_gases) + len(self.forwardmodel.atmosphere.inactive_gases)
    #     for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
    #         sumlog +=  log(self.forwardmodel.atmosphere.absorbing_gases_X[idx])
    #     for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
    #         sumlog += log(self.forwardmodel.atmosphere.inactive_gases_X[idx])
    #     mean = exp((1./n)*sumlog)
    #
    #     # create array of clr(absorbing_X + inactive_X)
    #     clr = []
    #     for idx, gasname in enumerate(self.forwardmodel.atmosphere.absorbing_gases):
    #         clr.append(log(self.forwardmodel.atmosphere.absorbing_gases_X[idx]/mean))
    #     for idx, gasname in enumerate(self.forwardmodel.atmosphere.inactive_gases):
    #         clr.append(log(self.forwardmodel.atmosphere.inactive_gases_X[idx]/mean))
    #
    #     return clr
    #
    # def get_mixing_ratio_inv_clr(self, clr):
    #
    #     # Convert clr(X) to X (centerd-log-ratio inverse transformation)
    #
    #     # transform back to simplex space
    #     clr_inv_tmp = []
    #     for i in range(self.forwardmodel.atmosphere.nallgases):
    #         clr_inv_tmp.append(exp(clr[i]))
    #     clr_inv_tmp = asarray(clr_inv_tmp)/sum(clr_inv_tmp) # closure
    #     clr_inv = clr_inv_tmp
    #
    #
    #     return clr_inv
