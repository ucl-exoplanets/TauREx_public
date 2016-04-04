'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Parameters class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

from ConfigParser import SafeConfigParser
import numpy as np
import logging
import os
import inspect
import subprocess

try:
    from mpi4py import MPI
    MPIimport = True
except ImportError:
    MPIimport = False
    pass
# MPI support
if MPIimport:
    MPIrank     = MPI.COMM_WORLD.Get_rank()
    MPIsize     = MPI.COMM_WORLD.Get_size()
else:
    MPIrank     = 0
    MPIsize     = 0

from library_constants import *

class parameters(object):

    def __init__(self, parfile='Parfiles/default.par'):
        '''

        a parameter file is parsed and initial parameter values are set.
        to add a new parameter edit this file and the input .par file.
        '''

        #config file parser
        if parfile:
            self.parser = SafeConfigParser()
            self.parser.readfp(open(parfile, 'rb'))

        self.parfile = parfile
        self.default_parser = SafeConfigParser()
        self.default_parser.read('Parfiles/default.par')
        self.default_parser.sections()
        self.verbose = self.getpar('General', 'verbose', 'bool')
        self.verbose_all_threads = self.getpar('General', 'verbose_all_threads', 'bool')

        if len(logging.getLogger().handlers) == 0: # be sure to load only one logging handler
            # configure logging instance
            if MPIrank == 0:
                if not os.path.isdir(self.getpar('Output','path')):
                        logging.info('Create folder Output')
                        os.makedirs(self.getpar('Output','path'))
                logging.basicConfig(filename=os.path.join(self.getpar('Output','path'), 'taurex.log'),
                                    level=logging.DEBUG)

            if (MPIrank == 0 and not self.verbose_all_threads) or self.verbose_all_threads:

                # define a Handler which writes INFO messages or higher to the sys.stderr
                self.console = logging.StreamHandler()
                self.console.setLevel(logging.DEBUG)
                formatter = logging.Formatter('%(asctime)s - Thread ' + str(MPIrank) + ' - %(levelname)s - %(message)s')
                self.console.setFormatter(formatter)
                logging.getLogger().addHandler(self.console)
                logging.info('Log started. Verbose for all threads: %s' % self.verbose_all_threads)

        self.version = subprocess.check_output(["git", "describe"])
        logging.info('Running TauREx %s' % self.version)

        logging.info('Initialise parameters object')

        # list of all molecules for which we have cross sections
        self.all_absorbing_gases = ['H2O', 'HCN', 'CH4', 'CO2', 'CO', 'NH3', 'C2H2']
        # list of all inactive gases we take care of
        self.all_inactive_gases = ['He', 'H2', 'N2']

        # section General
        self.gen_manual_waverange  = self.getpar('General','manual_waverange', 'bool')
        self.gen_wavemin           = self.getpar('General','wavemin', 'float')
        self.gen_wavemax           = self.getpar('General','wavemax', 'float')
        self.gen_type              = self.getpar('General','type')
        self.gen_ace              = self.getpar('General','ace', 'bool')
        self.gen_compile_cpp       = self.getpar('General','compile_cpp', 'bool')
        self.gen_run_gui           = False

        # section Input
        self.in_spectrum_file      = self.getpar('Input','spectrum_file')
        self.in_spectrum_db        = self.getpar('Input', 'spectrum_db')
        self.in_use_ATMfile        = self.getpar('Input','use_ATMfile', 'bool')
        self.in_atm_file           = self.getpar('Input','atm_file')
        self.in_xsec_alltemp       = False
        self.in_xsec_path          = self.getpar('Input','xsec_path')
        self.in_xsec_dnu           = self.getpar('Input','xsec_dnu', 'float')
        self.in_cia_path           = self.getpar('Input','cia_path')
        self.in_star_path          = self.getpar('Input','star_path')

        # section Output
        self.out_path              = self.getpar('Output','path')
        self.out_save_plots        = self.getpar('Output','save_plots', 'bool')
        self.out_sigma_spectrum        = self.getpar('Output', 'sigma_spectrum', 'bool')
        self.out_sigma_spectrum_frac   = self.getpar('Output', 'sigma_spectrum_frac', 'float')

        # section Star
        self.star_radius           = self.getpar('Star', 'radius', 'float')    *RSOL
        self.star_temp             = self.getpar('Star','temp', 'float')

        # section Planet
        self.planet_class          = self.getpar('Planet','class')
        self.planet_radius         = self.getpar('Planet', 'radius', 'float')*RJUP
        self.planet_mass           = self.getpar('Planet', 'mass', 'float')*MJUP

        # section Atmosphere
        self.atm_nlayers            = self.getpar('Atmosphere', 'nlayers', 'int')
        self.atm_max_pres           = self.getpar('Atmosphere', 'max_pressure', 'float')
        self.atm_min_pres           = self.getpar('Atmosphere', 'min_pressure', 'float')

        self.atm_tp_type              = self.getpar('Atmosphere', 'tp_type')
        self.atm_tp_iso_temp          = self.getpar('Atmosphere', 'tp_iso_temp', 'float' )
        self.atm_tp_guillot_T_irr     = self.getpar('Atmosphere', 'tp_guillot_T_irr', 'float' )
        self.atm_tp_guillot_kappa_irr = self.getpar('Atmosphere', 'tp_guillot_kappa_irr', 'float' )
        self.atm_tp_guillot_kappa_v1  = self.getpar('Atmosphere', 'tp_guillot_kappa_v1', 'float' )
        self.atm_tp_guillot_kappa_v2  = self.getpar('Atmosphere', 'tp_guillot_kappa_v2', 'float' )
        self.atm_tp_guillot_alpha     = self.getpar('Atmosphere', 'tp_guillot_alpha', 'float' )
        self.atm_tp_2point_T_surf     = self.getpar('Atmosphere', 'tp_2point_T_surf', 'float')
        self.atm_tp_2point_T_trop_diff= self.getpar('Atmosphere', 'tp_2point_T_trop_diff', 'float')
        self.atm_tp_2point_P_trop     = self.getpar('Atmosphere', 'tp_2point_P_trop', 'float')
        self.atm_tp_corr_length       = self.getpar('Atmosphere', 'tp_corr_length','float')
        

        self.atm_active_gases       = [gas.upper() for gas in self.getpar('Atmosphere','active_gases', 'list-str')]
        self.atm_active_gases_mixratios = self.getpar('Atmosphere','active_gases_mixratios', 'list-float')
        self.atm_inactive_gases     = [gas.upper() for gas in self.getpar('Atmosphere','inactive_gases', 'list-str')]
        self.atm_inactive_gases_mixratios = self.getpar('Atmosphere','inactive_gases_mixratios', 'list-float')

        self.atm_mu                 = self.getpar('Atmosphere', 'mu', 'float')*AMU
        self.atm_couple_mu          = self.getpar('Atmosphere', 'couple_mu', 'bool')

        self.atm_rayleigh           = self.getpar('Atmosphere','rayleigh', 'bool')
        self.atm_cia                = self.getpar('Atmosphere','cia', 'bool')
        self.atm_cia_pairs          = [pair.upper() for pair in self.getpar('Atmosphere','cia_pairs', 'list-str')]
        self.atm_clouds             = self.getpar('Atmosphere','clouds', 'bool')
        self.atm_cld_topP           = self.getpar('Atmosphere','cld_topP', 'float')
        self.atm_ace_metallicity    = self.getpar('Atmosphere', 'ace_metallicity', 'float')
        self.atm_ace_co             = self.getpar('Atmosphere', 'ace_co', 'float')

        # section Venot
        self.ven_load = self.getpar('Venot', 'load', 'bool')
        self.ven_TP_profile_path = self.getpar('Venot', 'TP_profile_path')
        self.ven_mol_profile_path = self.getpar('Venot', 'mol_profile_path')
        self.ven_exclude_mol = [mol.upper() for mol in self.getpar('Venot','exclude_mol', 'list-str')]

        # Section Fit

        self.fit_transmission      = self.getpar('Fitting','transmission', 'bool')
        self.fit_emission          = self.getpar('Fitting', 'emission', 'bool')
        self.fit_emission_stage2   = self.getpar('Fitting', 'emission_stage2', 'bool')

        # misc
        self.fit_couple_mu           = self.getpar('Fitting','couple_mu', 'bool')
        self.fit_inactive_mu_rescale = self.getpar('Fitting','inactive_mu_rescale', 'bool')
        self.fit_X_log               = self.getpar('Fitting','X_log', 'bool')
        self.fit_clr_trans           = self.getpar('Fitting','clr_trans', 'bool')

        # fit / fix parameters
        self.fit_fit_active          = self.getpar('Fitting', 'fit_active', 'bool')
        self.fit_fit_inactive        = self.getpar('Fitting', 'fit_inactive', 'bool')
        self.fit_fit_temp            = self.getpar('Fitting', 'fit_temp', 'bool')
        self.fit_fit_mu              = self.getpar('Fitting', 'fit_mu', 'bool')
        self.fit_fit_radius          = self.getpar('Fitting', 'fit_radius', 'bool')
        self.fit_fit_P0              = self.getpar('Fitting', 'fit_P0', 'bool')
        self.fit_fit_clouds_topP     = self.getpar('Fitting', 'fit_clouds_topP', 'bool')
        self.fit_fit_ace_metallicity = self.getpar('Fitting', 'fit_ace_metallicity', 'bool')
        self.fit_fit_ace_co          = self.getpar('Fitting', 'fit_ace_co', 'bool')

        # prior bounds
        self.fit_X_active_bounds        = self.getpar('Fitting', 'X_active_bounds', 'list-float')
        self.fit_X_inactive_bounds      = self.getpar('Fitting', 'X_inactive_bounds', 'list-float')
        self.fit_clr_bounds             = self.getpar('Fitting', 'clr_bounds', 'list-float')
        self.fit_mu_bounds              = self.getpar('Fitting', 'mu_bounds', 'list-float')
        self.fit_radius_bounds          = self.getpar('Fitting', 'radius_bounds', 'list-float')
        self.fit_P0_bounds              = self.getpar('Fitting', 'P0_bounds', 'list-float')
        self.fit_clouds_topP_bounds     = self.getpar('Fitting', 'clouds_topP_bounds', 'list-float')
        self.fit_ace_metallicity_bounds = self.getpar('Fitting', 'ace_metallicity_bounds', 'list-float')
        self.fit_ace_co_bounds          = self.getpar('Fitting', 'ace_co_bounds', 'list-float')

        self.fit_tp_iso_bounds               = self.getpar('Fitting', 'tp_iso_bounds', 'list-float')
        self.fit_tp_guillot_T_irr_bounds     = self.getpar('Fitting', 'tp_guillot_T_irr_bounds', 'list-float')
        self.fit_tp_guillot_kappa_irr_bounds = self.getpar('Fitting', 'tp_guillot_kappa_irr_bounds', 'list-float')
        self.fit_tp_guillot_kappa_v1_bounds  = self.getpar('Fitting', 'tp_guillot_kappa_v1_bounds', 'list-float')
        self.fit_tp_guillot_kappa_v2_bounds  = self.getpar('Fitting', 'tp_guillot_kappa_v2_bounds', 'list-float')
        self.fit_tp_guillot_alpha_bounds     = self.getpar('Fitting', 'tp_guillot_alpha_bounds', 'list-float')
        self.fit_hybrid_alpha_bounds         = self.getpar('Fitting', 'hybrid_alpha_bounds', 'list-float')

        # section Downhill
        self.downhill_run          = self.getpar('Downhill','run', 'bool')
        self.downhill_type         = self.getpar('Downhill', 'type')
        self.downhill_out_filename = self.getpar('Downhill','out_filename')

        # section MCMC
        self.mcmc_run              = self.getpar('MCMC','run', 'bool')
        self.mcmc_update_std       = self.getpar('MCMC','update_std', 'bool')
        self.mcmc_iter             = self.getpar('MCMC', 'iter', 'float')
        self.mcmc_burn             = self.getpar('MCMC','burn', 'float')
        self.mcmc_thin             = self.getpar('MCMC', 'thin', 'float')
        self.mcmc_verbose          = self.getpar('MCMC', 'verbose', 'bool')
        self.mcmc_progressbar      = self.getpar('MCMC', 'progressbar', 'bool')
        self.mcmc_out_filename     = self.getpar('MCMC','out_filename')

        # section Nest
        self.nest_run              = self.getpar('MultiNest','run', 'bool')
        self.nest_resume           = self.getpar('MultiNest','resume', 'bool')
        self.nest_verbose          = self.getpar('MultiNest','verbose', 'bool')
        self.nest_path             = self.getpar('MultiNest','nest_path') # @todo not used?
        self.nest_samp_eff         = self.getpar('MultiNest','sampling_eff')
        self.nest_nlive            = self.getpar('MultiNest','n_live_points', 'int')
        self.nest_max_iter         = self.getpar('MultiNest','max_iter', 'int')
        self.nest_multimodes       = self.getpar('MultiNest','multimodes')
        self.nest_max_modes        = self.getpar('MultiNest','max_modes', 'int')
        self.nest_const_eff        = self.getpar('MultiNest','const_eff', 'bool')
        self.nest_ev_tol           = self.getpar('MultiNest','evidence_tolerance','float')
        self.nest_mode_tol         = self.getpar('MultiNest', 'mode_tolerance', 'float')
        self.nest_imp_sampling     = self.getpar('MultiNest','imp_sampling', 'bool')
        self.nest_out_filename     = self.getpar('MultiNest','out_filename')

        #checking that either emission or transmisison is run
        if self.fit_emission and self.fit_transmission:
            logging.error('Transmission and emission cannot currently be run simultaneously')
            logging.error('change the fit_emission, fit_transmission parameters.' )

    def getpar(self, sec, par, type=None):

        # get parameter from user defined parser. If parameter is not found there, load the default parameter
        # the default parameter file parser is self.default_parser, defined in init

        try:

            if type == None:
                try:
                    return self.parser.get(sec, par)
                except:
                    return self.default_parser.get(sec, par)
            elif type == 'float':
                try:
                    return self.parser.getfloat(sec, par)
                except Exception,e:
                    return self.default_parser.getfloat(sec, par)

            elif type == 'bool':
                try:
                    return self.parser.getboolean(sec, par)
                except:
                    return self.default_parser.getboolean(sec, par)
            elif type == 'int':
                try:
                    return self.parser.getint(sec, par)
                except:
                    return self.default_parser.getint(sec, par)
            elif type == 'list-str':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [str(m).strip() for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [str(m).strip() for m in l]
            elif type == 'list-float':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [float(m) for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [float(m) for m in l]
            elif type == 'list-int':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [int(m) for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [int(m) for m in l]
            else:
                logging.error('Cannot set parameter %s in section %s. Parameter type %s not recognized. Set to None' (par, sec, type))
                return None
        except:
            logging.error('Cannot set parameter %s in section %s. Set to None' % (par, sec))
            return None

    def params_to_dict(self):

        # covert param variables to dictionary
        pr = {}
        for name in dir(self):
            value = getattr(self, name)
            if not name.startswith('__') and not inspect.ismethod(value) and \
                            name <> 'parser' and name <> 'default_parser' and name <> 'console':
                pr[name] = value
        return pr

