################################################
#class parameters 
#Parse parameter file eg. 'exonest.par' and initialise  
#parameters for run
################################################
from ConfigParser import SafeConfigParser
from numpy import genfromtxt,arange,size
from StringIO import StringIO

class parameters(object):
#instantiation
    def __init__(self, parfile):
        '''
        a parameter file is parsed and initial parameter values are set.  
        to add a new parameter edit this file and the input .par file.
        V1.0  - Definition - C. MacTavish, Apr 2012
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
        
        self.in_spectrum_file      = parser.get('Input','spectrum_file')
        self.in_use_ATMfile        = parser.getboolean('Input','use_ATMfile')
        self.in_atm_file           = parser.get('Input','atm_file')
        self.in_abs_path           = parser.get('Input','abs_path')
        self.in_abs_files          = parser.get('Input','abs_files')

        self.in_include_rad        = parser.getboolean('Input','include_rad')
        self.in_rad_file           = parser.get('Input','rad_file')
        self.in_include_cia        = parser.getboolean('Input','include_cia')
        self.in_cia_file           = parser.get('Input','cia_file')  
        self.in_include_cld        = parser.getboolean('Input','include_cld')
        self.in_cld_file           = parser.get('Input','cld_file')

        self.out_path              = parser.get('Output','path')
        self.out_file_prefix       = parser.get('Output','file_prefix')

        self.pre_run               = parser.getboolean('Preselector','run_pre')
        self.pre_cross_path        = parser.get('Preselector','cross_path')
        self.pre_speclib_path      = parser.get('Preselector','speclib_path')
        self.pre_pca_path          = parser.get('Preselector','pca_path')
        self.pre_conver2microns    = parser.getboolean('Preselector','convert2microns')
        self.pre_gen_speclib       = parser.getboolean('Preselector','generate_speclib')
        self.pre_mixing_ratios     = genfromtxt(StringIO(parser.get('Preselector', 'mixing_ratio')), delimiter = ',')
        self.pre_gen_pca           = parser.getboolean('Preselector','generate_pca')
        
        self.star_radius           = parser.getfloat('Star', 'radius')    *RSOL
        self.star_temp             = parser.getfloat('Star','temp')
        
        self.planet_radius         = parser.getfloat('Planet', 'radius')  *RJUP
        self.planet_sma            = parser.getfloat('Planet', 'sma')     *AU
        self.planet_grav           = parser.getfloat('Planet', 'grav')
        self.planet_albedo         = parser.getfloat('Planet','albedo')
        self.planet_temp           = parser.getfloat('Planet', 'temp')
        self.planet_mu             = parser.getfloat('Planet', 'mu')      *AMU
        self.planet_molec          = parser.get('Planet','molec')
        
        self.tp_var_atm            = parser.getboolean('T-P profile','var_atm')
        self.tp_num_scale          = parser.getint('T-P profile', 'num_scaleheights')
        self.tp_atm_levels         = parser.getint('T-P profile', 'atm_levels')
        self.tp_num_gas            = parser.getint('T-P profile', 'num_gas')
        self.tp_var_temp           = parser.getboolean('T-P profile', 'var_temp')
        self.tp_var_pres           = parser.getboolean('T-P profile', 'var_pres')
        self.tp_max_pres           = parser.getfloat('T-P profile', 'atm_max_pressure')
        self.tp_var_mix            = parser.getboolean('T-P profile', 'var_mix')

        self.fit_spec_res          = parser.getfloat('Fitting','spec_res')
        self.fit_transmission      = parser.getboolean('Fitting','transmission')
        self.fit_emission          = parser.getboolean('Fitting', 'emission')
        self.fit_param_free        = genfromtxt(StringIO(parser.get('Fitting', 'param_free')), delimiter = ',')
        self.fit_param_free_T      = arange(self.fit_param_free[0],dtype=int)
        self.fit_param_free_X      = arange(self.fit_param_free[0],self.fit_param_free[1]+self.fit_param_free[0],dtype=int)
        self.fit_T_up              = parser.getfloat('Fitting','T_up')
        self.fit_T_low             = parser.getfloat('Fitting','T_low')
        self.fit_X_up              = parser.getfloat('Fitting','X_up')
        self.fit_X_low             = parser.getfloat('Fitting','X_low')
        
        self.mcmc_update_std       = parser.getboolean('MCMC','update_std')
        self.mcmc_iter             = parser.getfloat('MCMC', 'iter')
        self.mcmc_burn             = parser.getfloat('MCMC','burn')
        self.mcmc_thin             = parser.getfloat('MCMC', 'thin')

        self.nest_resume           = parser.getboolean('MultiNest','resume')
        self.nest_verbose          = parser.getboolean('MultiNest','verbose')
        self.nest_samp_eff         = parser.get('MultiNest','sampling_eff')
        self.nest_nlive            = parser.getint('MultiNest','n_live_points')
        self.nest_max_iter         = parser.getint('MultiNest','max_iter')
        self.nest_imp_sampling     = parser.getboolean('MultiNest','imp_sampling')
        
        
        
    def list(self,name=None):
        if name is None:
            return dir(self)[2:-1]
        else:
            lst = dir(self)
            return filter(lambda k: name in k, lst)
        
    def __getattribute__(self,name):
        return object.__getattribute__(self, name)
    
    def __getitem__(self,name):
        return self.__dict__[name]

