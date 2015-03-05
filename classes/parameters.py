################################################
#class parameters 
#Parse parameter file eg. 'exonest.par' and initialise  
#parameters for run
################################################
from base import base
from ConfigParser import SafeConfigParser
import numpy as np
from numpy import genfromtxt,arange,size
from StringIO import StringIO
import ast
import logging
import os

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


#conversion constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg


class parameters(base):
#instantiation
    def __init__(self, parfile='Parfiles/default.par'):
        '''
        a parameter file is parsed and initial parameter values are set.  
        to add a new parameter edit this file and the input .par file.
        V1.0  - Definition - I. Waldmann, Apr 2013
        '''

        #config file parser
        if parfile:
            self.parser = SafeConfigParser()
            self.parser.readfp(open(parfile, 'rb'))

        self.default_parser = SafeConfigParser()
        self.default_parser.read('Parfiles/default.par')
        self.default_parser.sections()
        self.verbose = self.getpar('General', 'verbose', 'bool')
        self.verbose_all_threads = self.getpar('General', 'verbose_all_threads', 'bool')

        if len(logging.getLogger().handlers) == 0: # be sure to load only one logging handler
            # configure logging instance
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
        self.all_absorbing_gases = ['1H2-16O', '1H-12C-14N', '12C-1H4', '12C-16O2', '12C-16O', '14N-1H3',
                                    '28Si-16O', '48Ti-16O', '14N-16O'] # VO removed due to lowest T = 1000 in cross sections

        # list of all inactive gases we take care of
        self.all_inactive_gases = ['He', 'H2', 'N2']

        self.trans_cpp             = self.getpar('General', 'trans_cpp', 'bool')
        
        self.gen_manual_waverange  = self.getpar('General','manual_waverange', 'bool')
        self.gen_wavemin           = self.getpar('General','wavemin', 'float')
        self.gen_wavemax           = self.getpar('General','wavemax', 'float')
        self.gen_spec_res          = self.getpar('General','spec_res', 'float')
        self.gen_type              = self.getpar('General','type')
        self.gen_run_gui           = False

        self.in_spectrum_file      = self.getpar('Input','spectrum_file')
        if self.in_spectrum_file == 'False':
            self.in_spectrum_file = False
        self.in_use_ATMfile        = self.getpar('Input','use_ATMfile', 'bool')
        self.in_atm_file           = self.getpar('Input','atm_file')
        self.in_abs_path           = self.getpar('Input','abs_path')
        self.in_convert2microns    = self.getpar('Input','convert2microns', 'bool')
        self.in_abs_files          = self.getpar('Input','__legacy__abs_files')
        
        self.in_tempres            = self.getpar('Input','tempres', 'float')
        self.in_star_path          = self.getpar('Input','star_path')
        self.in_include_rad        = self.getpar('Input','include_rad', 'bool')
        self.in_rad_file           = self.getpar('Input','rad_file')
        self.in_include_cia        = self.getpar('Input','include_cia', 'bool')
        self.in_cia_file           = self.getpar('Input','cia_file')
        self.in_create_sigma_rayleigh = self.getpar('Input','create_sigma_rayleigh', 'bool')

        self.out_path              = self.getpar('Output','path')
        self.out_file_prefix       = self.getpar('Output','file_prefix')
        self.out_dump_internal     = self.getpar('Output','dump_internal', 'bool')
        self.out_internal_name     = self.getpar('Output','internal_name')
        self.out_save_plots        = self.getpar('Output','save_plots', 'bool')
        self.out_plot_contour        = self.getpar('Output','plot_contour', 'bool')

        self.star_radius           = self.getpar('Star', 'radius', 'float')    *RSOL
        self.star_temp             = self.getpar('Star','temp', 'float')
        
        self.planet_name           = self.getpar('Planet', 'name')
        self.planet_radius         = self.getpar('Planet', 'radius', 'float')  *RJUP
        self.planet_mass           = self.getpar('Planet', 'mass', 'float')     *MJUP
        self.planet_sma            = self.getpar('Planet', 'sma', 'float')     *AU
        self.planet_albedo         = self.getpar('Planet','albedo', 'float')
        self.planet_temp           = self.getpar('Planet', 'temp', 'float')
        self.planet_mu             = self.getpar('Planet', 'mu', 'float')      *AMU
        self.planet_molec          = self.getpar('Planet','molecules', 'list-str')  #genfromtxt(StringIO(self.getpar('Planet','molecules')),delimiter = ',',dtype='str',autostrip=True)
        self.planet_mixing         = self.getpar('Planet','mixing_ratios', 'list-float') #genfromtxt(StringIO(self.getpar('Planet','mixing_ratios')),delimiter = ',',dtype='str',autostrip=True)
        self.planet_inactive_gases     = self.getpar('Planet', 'inactive_gases', 'list-str')
        self.planet_inactive_gases_X   = self.getpar('Planet', 'inactive_gases_X', 'list-float')
        self.in_include_Rayleigh       = self.getpar('Planet','include_Rayleigh', 'bool')
        try:
            self.in_include_cld        = self.getpar('Planet','include_cld', 'bool')
            self.in_cld_params         = self.getpar('Planet','cld_params', 'list-float') #genfromtxt(StringIO(self.getpar('Planet','cld_params')),delimiter = ',',dtype='str',autostrip=True)
            self.in_cld_m              = np.float(self.in_cld_params[0])
            self.in_cld_a              = np.float(self.in_cld_params[1])
            self.in_cld_pressure       = self.getpar('Planet','cld_pressure', 'list-float') #np.float(genfromtxt(StringIO(self.getpar('Planet','cld_pressure')),delimiter = ',',dtype='str',autostrip=True))
            self.in_cld_file           = self.getpar('Planet','cld_file')
        except:
            self.in_include_cld        = False
            pass




        self.tp_var_atm            = self.getpar('T-P profile','var_atm', 'bool')
        self.tp_num_scale          = self.getpar('T-P profile', 'num_scaleheights', 'int')
        self.tp_atm_levels         = self.getpar('T-P profile', 'atm_levels', 'int')
        #self.tp_num_gas           = self.getpar('T-P profile', 'num_gas', 'int') # deprecated
        self.tp_var_temp           = self.getpar('T-P profile', 'var_temp', 'bool')
        self.tp_var_pres           = self.getpar('T-P profile', 'var_pres', 'bool')
        self.tp_max_pres           = self.getpar('T-P profile', 'atm_max_pressure', 'float')
        self.tp_var_mix            = self.getpar('T-P profile', 'var_mix', 'bool')
        self.tp_type               = self.getpar('T-P profile', 'profile_type')


        try:
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
        except:
            self.pre_run               = False
            pass

        try:
            self.fit_transmission      = self.getpar('Fitting','transmission', 'bool')
            self.fit_emission          = self.getpar('Fitting', 'emission', 'bool')
            self.fit_emission_stage2   = self.getpar('Fitting', 'emission_stage2', 'bool')
            self.fit_param_free        = self.getpar('Fitting', 'param_free', 'list-float') # currently not used
            self.fit_param_free_T      = arange(self.fit_param_free[0],dtype=int)
            self.fit_param_free_X      = arange(self.fit_param_free[0],self.fit_param_free[1]+self.fit_param_free[0],dtype=int)
            self.fit_fix_temp          = self.getpar('Fitting', 'fix_temp', 'bool')
            self.fit_fix_P0            = self.getpar('Fitting', 'fix_P0', 'bool')
            self.fit_fix_radius        = self.getpar('Fitting', 'fix_radius', 'bool')
            self.fit_fix_mu            = self.getpar('Fitting', 'fix_mu', 'bool')
            self.fit_fix_inactive      = self.getpar('Fitting', 'fix_inactive', 'bool')
            self.fit_couple_mu         = self.getpar('Fitting', 'couple_mu', 'bool')
            self.fit_T_up              = self.getpar('Fitting','T_up', 'float')
            self.fit_T_low             = self.getpar('Fitting','T_low', 'float')
            self.fit_radius_up         = self.getpar('Fitting','radius_up', 'float')  # in RJUP
            self.fit_radius_low        = self.getpar('Fitting','radius_low', 'float') # in RJUP
            self.fit_P0_up             = self.getpar('Fitting','P0_up', 'float')  # in Pascal
            self.fit_P0_low            = self.getpar('Fitting','P0_low', 'float') # in Pascal
            self.fit_mu_up             = self.getpar('Fitting','mu_up', 'float') # in AMU
            self.fit_mu_low            = self.getpar('Fitting','mu_low', 'float') # in AMU
            self.fit_X_up              = self.getpar('Fitting','X_up', 'float')
            self.fit_X_low             = self.getpar('Fitting','X_low', 'float')
            self.fit_X_inactive_up     = self.getpar('Fitting','X_inactive_up', 'float')
            self.fit_X_inactive_low    = self.getpar('Fitting','X_inactive_low', 'float')
            self.fit_clr_trans    = self.getpar('Fitting','clr_trans', 'bool')
        except:
            self.fit_transmission      = False
            self.fit_emission          = False
            pass
        try: 
            self.downhill_run          = self.getpar('Downhill','run', 'bool')
            self.downhill_type         = self.getpar('Downhill', 'type')
#         self.downhill_options      = ast.literal_eval(str(self.getpar('Downhill','options')))
        except:
            self.downhill_run          = False
            pass
        
        try:
            self.mcmc_run              = self.getpar('MCMC','run', 'bool')
            self.mcmc_update_std       = self.getpar('MCMC','update_std', 'bool')
            self.mcmc_iter             = self.getpar('MCMC', 'iter', 'float')
            self.mcmc_burn             = self.getpar('MCMC','burn', 'float')
            self.mcmc_thin             = self.getpar('MCMC', 'thin', 'float')
            self.mcmc_verbose          = self.getpar('MCMC', 'verbose', 'bool')
            self.mcmc_progressbar          = self.getpar('MCMC', 'progressbar', 'bool')
        except:
            self.mcmc_run              = False
            pass

        try:
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
            self.nest_imp_sampling     = self.getpar('MultiNest','imp_sampling', 'bool')
            self.nest_cluster_analysis = self.getpar('MultiNest','cluster_analysis', 'bool')
        except:
            self.nest_run              = False
            pass

        try:
            self.grid_res              = self.getpar('Grid','resolution', 'list-float')
            self.grid_snr              = self.getpar('Grid','snr', 'list-float')
        except:
            pass

        try:
            self.clean_run             = self.getpar('Housekeeping','run','bool')
            self.clean_script          = self.getpar('Housekeeping','script_name')
            self.clean_save_used_params= self.getpar('Housekeeping','save_used_params','bool')
        except:
            self.clean_run             = False
            pass


        #####################################################################
        #additional parser commands. Checking parameter compatibility and stuff 
        
        #checking that either emission or transmisison is run
        if self.fit_emission and self.fit_transmission:
            print 'Error: transmission and emission cannot currently be run simultaneously'
            print 'change the fit_emission, fit_transmission parameters.'  
        if self.fit_emission: self.fit_transmission = False
        if self.fit_transmission: self.fit_emission = False

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