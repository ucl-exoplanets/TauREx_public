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

class parameters(base):
#instantiation
    def __init__(self, parfile):
        '''
        a parameter file is parsed and initial parameter values are set.  
        to add a new parameter edit this file and the input .par file.
        V1.0  - Definition - I. Waldmann, Apr 2013
        '''
        
        #conversion constants
        RSOL  = 6.955e8         #stellar radius to m
        RJUP  = 6.9911e7        #jupiter radius to m
        REARTH= 6.371e3         #earth radius to m
        AU    = 1.49e11         #semi-major axis (AU) to m
        AMU   = 1.660538921e-27 #atomic mass to kg
        
        #config file parser
        parser = SafeConfigParser()
        parser.read(parfile)
        
        self.verbose               = parser.getboolean('General', 'verbose')
        self.trans_cpp             = parser.getboolean('General', 'trans_cpp')
        
        self.gen_manual_waverange  = parser.getboolean('General','manual_waverange')
        self.gen_wavemin           = parser.getfloat('General','wavemin')
        self.gen_wavemax           = parser.getfloat('General','wavemax')
        self.gen_spec_res          = parser.getfloat('General','spec_res')
        self.gen_type              = parser.get('General','type')
        
        self.in_spectrum_file      = parser.get('Input','spectrum_file')
        self.in_use_ATMfile        = parser.getboolean('Input','use_ATMfile')
        self.in_atm_file           = parser.get('Input','atm_file')
        self.in_abs_path           = parser.get('Input','abs_path')
        self.in_abs_files          = parser.get('Input','__legacy__abs_files')
        
        self.in_tempres            = parser.getfloat('Input','tempres')

        self.in_star_path          = parser.get('Input','star_path')
        
        self.in_include_rad        = parser.getboolean('Input','include_rad')
        self.in_rad_file           = parser.get('Input','rad_file')
        self.in_include_cia        = parser.getboolean('Input','include_cia')
        self.in_cia_file           = parser.get('Input','cia_file')  

        self.out_path              = parser.get('Output','path')
        self.out_file_prefix       = parser.get('Output','file_prefix')
        self.out_dump_internal     = parser.getboolean('Output','dump_internal')
        self.out_internal_name     = parser.get('Output','internal_name')
        self.out_save_plots        = parser.getboolean('Output','save_plots')

        self.star_radius           = parser.getfloat('Star', 'radius')    *RSOL
        self.star_temp             = parser.getfloat('Star','temp')
        
        self.planet_radius         = parser.getfloat('Planet', 'radius')  *RJUP
        self.planet_sma            = parser.getfloat('Planet', 'sma')     *AU
        self.planet_grav           = parser.getfloat('Planet', 'grav')
        self.planet_albedo         = parser.getfloat('Planet','albedo')
        self.planet_temp           = parser.getfloat('Planet', 'temp')
        self.planet_mu             = parser.getfloat('Planet', 'mu')      *AMU
        self.planet_molec          = genfromtxt(StringIO(parser.get('Planet','molecules')),delimiter = ',',dtype='str',autostrip=True)
        self.planet_mixing         = genfromtxt(StringIO(parser.get('Planet','mixing_ratios')),delimiter = ',',dtype='str',autostrip=True)
        
        try:
            self.in_include_cld        = parser.getboolean('Planet','include_cld')
            self.in_cld_params         = genfromtxt(StringIO(parser.get('Planet','cld_params')),delimiter = ',',dtype='str',autostrip=True)
            self.in_cld_m              = np.float(self.in_cld_params[0])
            self.in_cld_a              = np.float(self.in_cld_params[1])
            self.in_cld_pressure       = np.float(genfromtxt(StringIO(parser.get('Planet','cld_pressure')),delimiter = ',',dtype='str',autostrip=True))
            self.in_cld_file           = parser.get('Planet','cld_file')
        except:
            self.in_include_cld        = False
            pass
        
        self.tp_var_atm            = parser.getboolean('T-P profile','var_atm')
        self.tp_num_scale          = parser.getint('T-P profile', 'num_scaleheights')
        self.tp_atm_levels         = parser.getint('T-P profile', 'atm_levels')
        #self.tp_num_gas            = parser.getint('T-P profile', 'num_gas')
        self.tp_var_temp           = parser.getboolean('T-P profile', 'var_temp')
        self.tp_var_pres           = parser.getboolean('T-P profile', 'var_pres')
        self.tp_max_pres           = parser.getfloat('T-P profile', 'atm_max_pressure')
        self.tp_var_mix            = parser.getboolean('T-P profile', 'var_mix')


        try:
            self.pre_run               = parser.getboolean('Preselector','run_pre')
            self.pre_speclib_path      = parser.get('Preselector','speclib_path')
            self.pre_pca_path          = parser.get('Preselector','pca_path')
            self.pre_conver2microns    = parser.getboolean('Preselector','convert2microns')
            self.pre_gen_speclib       = parser.getboolean('Preselector','generate_speclib')
            self.pre_restrict_temp     = parser.getboolean('Preselector','restrict_temp')
            self.pre_temp_range        = genfromtxt(StringIO(parser.get('Preselector', 'temp_range')), delimiter = ',')
            self.pre_mixing_ratios     = genfromtxt(StringIO(parser.get('Preselector', 'mixing_ratio')), delimiter = ',')
            self.pre_gen_pca           = parser.getboolean('Preselector','generate_pca')
            self.pre_mask_thres        = parser.getfloat('Preselector','mask_thres')
            self.pre_mol_force_bool    = parser.getboolean('Preselector','mol_force_on')
            self.pre_mol_force         = genfromtxt(StringIO(parser.get('Preselector','mol_force')),delimiter = ',',dtype='str',autostrip=True)
        except:
            self.pre_run               = False
            pass

        try:
            self.fit_transmission      = parser.getboolean('Fitting','transmission')
            self.fit_emission          = parser.getboolean('Fitting', 'emission')
            self.fit_param_free        = genfromtxt(StringIO(parser.get('Fitting', 'param_free')), delimiter = ',')
            self.fit_param_free_T      = arange(self.fit_param_free[0],dtype=int)
            self.fit_param_free_X      = arange(self.fit_param_free[0],self.fit_param_free[1]+self.fit_param_free[0],dtype=int)
            self.fit_T_up              = parser.getfloat('Fitting','T_up')
            self.fit_T_low             = parser.getfloat('Fitting','T_low')
            self.fit_X_up              = parser.getfloat('Fitting','X_up')
            self.fit_X_low             = parser.getfloat('Fitting','X_low')
        except:
            self.fit_transmission      = False
            self.fit_emission          = False
            pass
        
        try:
            self.mcmc_run              = parser.getboolean('MCMC','run')
            self.mcmc_update_std       = parser.getboolean('MCMC','update_std')
            self.mcmc_iter             = parser.getfloat('MCMC', 'iter')
            self.mcmc_burn             = parser.getfloat('MCMC','burn')
            self.mcmc_thin             = parser.getfloat('MCMC', 'thin')
        except:
            self.mcmc_run              = False
            pass

        try:
            self.nest_run              = parser.getboolean('MultiNest','run')
            self.nest_resume           = parser.getboolean('MultiNest','resume')
            self.nest_verbose          = parser.getboolean('MultiNest','verbose')
            self.nest_path             = parser.get('MultiNest','nest_path')
            self.nest_samp_eff         = parser.get('MultiNest','sampling_eff')
            self.nest_nlive            = parser.getint('MultiNest','n_live_points')
            self.nest_max_iter         = parser.getint('MultiNest','max_iter')
            self.nest_imp_sampling     = parser.getboolean('MultiNest','imp_sampling')
        except:
            self.nest_run              = False
            pass
        
        
        #####################################################################
        #additional parser commands. Checking parameter compatibility and stuff 
        
        #checking that either emission or transmisison is run
        if self.fit_emission and self.fit_transmission:
            print 'Error: transmission and emission cannot currently be run simultaneously'
            print 'change the fit_emission, fit_transmission parameters.'  
        if self.fit_emission: self.fit_transmission = False
        if self.fit_transmission: self.fit_emission = False 
        
        

