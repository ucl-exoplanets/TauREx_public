'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Parameters class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

from base import base
from ConfigParser import SafeConfigParser
import numpy as np
from numpy import genfromtxt,arange,size
from StringIO import StringIO
import ast, logging, os


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

class parameters(base):
#instantiation
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
        self.gen_ACE              = self.getpar('General','ACE', 'bool')
        self.gen_compile_cpp       = self.getpar('General','compile_cpp', 'bool')
        self.gen_run_gui           = False

        # section Input
        self.in_spectrum_file      = self.getpar('Input','spectrum_file')
        if self.in_spectrum_file == 'False':
            self.in_spectrum_file = False
        self.in_use_ATMfile        = self.getpar('Input','use_ATMfile', 'bool')
        self.in_atm_file           = self.getpar('Input','atm_file')
        self.in_xsec_path          = self.getpar('Input','xsec_path')
        self.in_xsec_dnu           = self.getpar('Input','xsec_dnu', 'float')
        self.in_cia_path           = self.getpar('Input','cia_path')
        self.in_star_path          = self.getpar('Input','star_path')

        # section Output
        self.out_path              = self.getpar('Output','path')
        self.out_dump_internal     = self.getpar('Output','dump_internal', 'bool')
        self.out_internal_name     = self.getpar('Output','internal_name')
        self.out_save_plots        = self.getpar('Output','save_plots', 'bool')
        self.out_plot_contour      = self.getpar('Output','plot_contour', 'bool')
        self.out_plot_colour       = self.getpar('Output','plot_colour')

        # section Star
        self.star_radius           = self.getpar('Star', 'radius', 'float')    *RSOL
        self.star_temp             = self.getpar('Star','temp', 'float')

        # section Planet
        self.planet_class          = self.getpar('Planet','class')
        self.planet_radius         = self.getpar('Planet', 'radius', 'float')*RJUP
        self.planet_mass           = self.getpar('Planet', 'mass', 'float')*MJUP
        self.planet_temp           = self.getpar('Planet', 'temp', 'float')
        self.planet_mu             = self.getpar('Planet', 'mu', 'float')*AMU

        # section Atmosphere
        self.atm_num_scaleheights   = self.getpar('Atmosphere', 'num_scaleheights', 'int')
        self.atm_nlayers            = self.getpar('Atmosphere', 'nlayers', 'int')
        self.atm_max_pres           = self.getpar('Atmosphere', 'max_pressure', 'float')
        self.atm_min_pres           = self.getpar('Atmosphere', 'min_pressure', 'float')
        self.atm_tp_type            = self.getpar('Atmosphere', 'profile_type')
        self.atm_corrlength         = self.getpar('Atmosphere', 'corr_length','float')
        self.atm_active_gases       = [gas.upper() for gas in self.getpar('Atmosphere','active_gases', 'list-str')]
        self.atm_active_gases_mixratios = self.getpar('Atmosphere','active_gases_mixratios', 'list-float')
        self.atm_inactive_gases     = [gas.upper() for gas in self.getpar('Atmosphere','inactive_gases', 'list-str')]
        self.atm_inactive_gases_mixratios = self.getpar('Atmosphere','inactive_gases_mixratios', 'list-float')
        self.atm_couple_mu          = self.getpar('Atmosphere', 'couple_mu', 'bool')
        self.atm_rayleigh           = self.getpar('Atmosphere','rayleigh', 'bool')
        self.atm_cia                = self.getpar('Atmosphere','cia', 'bool')
        self.atm_cia_pairs          = [pair.upper() for pair in self.getpar('Atmosphere','cia_pairs', 'list-str')]
        self.atm_clouds             = self.getpar('Atmosphere','clouds', 'bool')
        self.atm_cld_params         = self.getpar('Atmosphere','cld_params', 'list-float')
        self.atm_cld_m              = self.atm_cld_params[0]
        self.atm_cld_a              = self.atm_cld_params[1]
        self.atm_cld_pressure       = self.getpar('Atmosphere','cld_pressure', 'list-float')
        self.atm_cld_lower_P        = self.atm_cld_pressure[0]
        self.atm_cld_upper_P        = self.atm_cld_pressure[1]
        self.atm_ace_He_abund_dex   = self.getpar('Atmosphere', 'ace_HE_abund_dex', 'float')
        self.atm_ace_C_abund_dex    = self.getpar('Atmosphere', 'ace_C_abund_dex', 'float')
        self.atm_ace_O_abund_dex    = self.getpar('Atmosphere', 'ace_O_abund_dex', 'float')
        self.atm_ace_N_abund_dex    = self.getpar('Atmosphere', 'ace_N_abund_dex', 'float')

        # section Venot
        self.ven_load = self.getpar('Venot', 'load', 'bool')
        self.ven_TP_profile_path = self.getpar('Venot', 'TP_profile_path')
        self.ven_mol_profile_path = self.getpar('Venot', 'mol_profile_path')
        self.ven_exclude_mol = [mol.upper() for mol in self.getpar('Venot','exclude_mol', 'list-str')]

        # section Preselecto/100r
        self.pre_run               = self.getpar('Preselector','run_pre', 'bool')
        self.pre_speclib_path      = self.getpar('Preselector','speclib_path')
        self.pre_pca_path          = self.getpar('Preselector','pca_path')
        self.pre_gen_speclib       = self.getpar('Preselector','generate_speclib', 'bool')
        self.pre_restrict_temp     = self.getpar('Preselector','restrict_temp', 'bool')
        self.pre_temp_range        = self.getpar('Preselector', 'temp_range', 'list-float')  #genfromtxt(StringIO(self.getpar('Preselector', 'temp_range')), delimiter = ',')
        self.pre_mixing_ratios     = self.getpar('Preselector', 'mixing_ratio', 'list-float')  #genfromtxt(StringIO(self.getpar('Preselector', 'mixing_ratio')), delimiter = ',')
        self.pre_gen_pca           = self.getpar('Preselector','generate_pca', 'bool')
        self.pre_mask_thres        = self.getpar('Preselector','mask_thres', 'float')
        self.pre_mol_force_bool    = self.getpar('Preselector','mol_force_on', 'bool')
        self.pre_mol_force         = self.getpar('Preselector', 'mol_force', 'list-str')  # genfromtxt(StringIO(self.getpar('Preselector','mol_force')),delimiter = ',',dtype='str',autostrip=True)

        # Section Fit
        self.fit_transmission      = self.getpar('Fitting','transmission', 'bool')
        self.fit_emission          = self.getpar('Fitting', 'emission', 'bool')
        self.fit_emission_stage2   = self.getpar('Fitting', 'emission_stage2', 'bool')

        self.fit_couple_mu    = self.getpar('Fitting','couple_mu', 'bool')
        self.fit_X_log        = self.getpar('Fitting','X_log', 'bool')
        self.fit_clr_trans    = self.getpar('Fitting','clr_trans', 'bool')

        self.fit_fit_active          = self.getpar('Fitting', 'fit_active', 'bool')
        self.fit_fit_inactive        = self.getpar('Fitting', 'fit_inactive', 'bool')
        self.fit_fit_temp            = self.getpar('Fitting', 'fit_temp', 'bool')
        self.fit_fit_mu              = self.getpar('Fitting', 'fit_mu', 'bool')
        self.fit_fit_radius          = self.getpar('Fitting', 'fit_radius', 'bool')
        self.fit_fit_P0              = self.getpar('Fitting', 'fit_P0', 'bool')
        self.fit_fit_clouds_lower_P  = self.getpar('Fitting', 'fit_clouds_lower_P', 'bool')
        self.fit_fit_clouds_upper_P  = self.getpar('Fitting', 'fit_clouds_upper_P', 'bool')
        self.fit_fit_clouds_m        = self.getpar('Fitting', 'fit_clouds_m', 'bool')
        self.fit_fit_clouds_a        = self.getpar('Fitting', 'fit_clouds_a', 'bool')
        self.fit_fit_He_abund_dex    = self.getpar('Fitting', 'fit_He_abund_dex', 'bool')
        self.fit_fit_C_abund_dex     = self.getpar('Fitting', 'fit_C_abund_dex', 'bool')
        self.fit_fit_O_abund_dex     = self.getpar('Fitting', 'fit_O_abund_dex', 'bool')
        self.fit_fit_N_abund_dex     = self.getpar('Fitting', 'fit_N_abund_dex', 'bool')

        self.fit_X_active_bounds       = self.getpar('Fitting', 'X_active_bounds', 'list-float')
        self.fit_X_inactive_bounds     = self.getpar('Fitting', 'X_inactive_bounds', 'list-float')
        self.fit_clr_bounds            = self.getpar('Fitting', 'clr_bounds', 'list-float')
        self.fit_T_bounds              = self.getpar('Fitting', 'T_bounds', 'list-float')
        self.fit_mu_bounds             = self.getpar('Fitting', 'mu_bounds', 'list-float')
        self.fit_radius_bounds         = self.getpar('Fitting', 'radius_bounds', 'list-float')
        self.fit_P0_bounds             = self.getpar('Fitting', 'P0_bounds', 'list-float')
        self.fit_clouds_lower_P_bounds = self.getpar('Fitting', 'clouds_lower_P_bounds', 'list-float')
        self.fit_clouds_upper_P_bounds = self.getpar('Fitting', 'clouds_upper_P_bounds', 'list-float')
        self.fit_clouds_a_bounds       = self.getpar('Fitting', 'clouds_a_bounds', 'list-float')
        self.fit_clouds_m_bounds       = self.getpar('Fitting', 'clouds_m_bounds', 'list-float')
        self.fit_ace_He_abund_dex_bounds = self.getpar('Fitting', 'ace_He_abund_dex_bounds', 'list-float')
        self.fit_ace_C_abund_dex_bounds  = self.getpar('Fitting', 'ace_C_abund_dex_bounds', 'list-float')
        self.fit_ace_O_abund_dex_bounds  = self.getpar('Fitting', 'ace_O_abund_dex_bounds', 'list-float')
        self.fit_ace_N_abund_dex_bounds  = self.getpar('Fitting', 'ace_N_abund_dex_bounds', 'list-float')

        self.fit_hybrid_alpha_l    = self.getpar('Fitting', 'hybrid_alpha_low', 'float')
        self.fit_hybrid_alpha_h    = self.getpar('Fitting', 'hybrid_alpha_high', 'float')

        # section Downhill
        self.downhill_run          = self.getpar('Downhill','run', 'bool')
        self.downhill_type         = self.getpar('Downhill', 'type')

        # section MCMC
        self.mcmc_run              = self.getpar('MCMC','run', 'bool')
        self.mcmc_update_std       = self.getpar('MCMC','update_std', 'bool')
        self.mcmc_iter             = self.getpar('MCMC', 'iter', 'float')
        self.mcmc_burn             = self.getpar('MCMC','burn', 'float')
        self.mcmc_thin             = self.getpar('MCMC', 'thin', 'float')
        self.mcmc_verbose          = self.getpar('MCMC', 'verbose', 'bool')
        self.mcmc_progressbar      = self.getpar('MCMC', 'progressbar', 'bool')

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
        self.nest_cluster_analysis = self.getpar('MultiNest','cluster_analysis', 'bool')

        # section Housekeeping ???
        try:
            self.clean_run             = self.getpar('Housekeeping','run','bool')
            self.clean_script          = self.getpar('Housekeeping','script_name')
            self.clean_save_used_params= self.getpar('Housekeeping','save_used_params','bool')
        except:
            self.clean_run             = False
            pass

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