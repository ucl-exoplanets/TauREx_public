################################################
#class profile
#
# Reads in all relevant data and performs pre-processing 
# such as sorting, grid-interpolations etc. 
#
# Input: -parameter object
#        -data object
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



class profile(base):


    def __init__(self,params,data):
        
        self.params = params
        
        #constants
        self.boltzmann = 1.3806488e-23# m2 kg s-2 K-1
        
        #derived values
        self.scaleheight = self.get_scaleheight(params.planet_temp, params.planet_grav, params.planet_mu)
     
        
        
        if params.tp_var_atm:
            self.nlayers = int(params.tp_atm_levels)
            self.ngas    = int(data.ngas)
            self.pta     = self.setup_pta_grid()
            self.P       = self.pta[:,0]
            self.T       = self.pta[:,1]
            self.Z       = self.pta[:,2]
            self.X       = np.zeros((self.ngas,self.nlayers))
            self.X      += 1e-5  #setting up initial mixing ratios
            self.rho     = self.get_rho()
            
        else:           
            self.nlayers = data.nlayers
            self.pta     = data.pta
            self.P       = self.pta[:,0]
            self.T       = data.pta[:,1]
            self.Z       = self.pta[:,2]
            self.X       = data.X       
            self.rho     = self.get_rho(T=self.T,P=self.P)
            
        T,P = self.TP_profile([1e3,10.0], [1800,1400,1600])
        
        pl.figure(200)
        pl.plot(T,P)
        pl.yscale('log')
        pl.gca().invert_yaxis()
        pl.show()
        exit()
            
                

    #class methods 
            
    def get_scaleheight(self,T_aver,surf_g,mmw):
        
        return (self.boltzmann*T_aver)/(mmw*surf_g)
        
        
    # @profile #line-by-line profiling decorator
    def setup_pta_grid(self):
        #calculate pressure, temperature, altitude grid
        
        MAX_P    = self.params.tp_max_pres
        N_SCALE  = self.params.tp_num_scale
        N_LAYERS = self.nlayers
        
        max_z    = N_SCALE * self.scaleheight
        
#         dz       = max_z / N_LAYERS
        
        #generatinng altitude-pressure array
        PTA_arr = np.zeros((N_LAYERS,3))
        PTA_arr[:,2] = np.linspace(0,max_z,num=N_LAYERS)
        PTA_arr[:,0] = MAX_P * np.exp(-PTA_arr[:,2]/self.scaleheight)
        PTA_arr[:,1] = self.params.planet_temp

        
        return PTA_arr
        
    def setup_prior_bounds(self):
        #partially to be moved to parameter file i guess
        
        pass

    def setup_parameter_grid(self, transmission=False, emission=False):
        #I AM SAMPSAN I AM SAMPSAN
        
        INDEX  = []
        PARAMS = []
        
        #setting up mixing ratios for individual gases
        cgas = 0
        for i in self.ngas:
            PARAMS.append(1e-5)
            cgas += 1
        
        INDEX.append(cgas)
        
        #setting up temperature parameters
        T_mean = (self.params.T_up - self.params.T_low)/2.0
        num_T_params = 3
        
        ctemp = 0
        if transmission:
            ctemp +=1
            PARAMS.append(self.params.planet_temp)
        if emission:
            for i in range(num_T_params):
                PARAMS.append(T_mean)
                ctemp += 1
            for i in range(num_T_params-1):
                
                
            #UNFINISHED!!
        
        

    # @profile #line-by-line profiling decorator
    def get_rho(self,T=None,P=None):
        #calculate atmospheric densities for given temperature and pressure
        
        if P is None:
            P = self.P
        if T is None:
            T = self.params.planet_temp
            
        return  (P)/(self.boltzmann*T)   
        
    # def cast_FIT_array(self,FIT,):
    
    def TP_profile(self,PARAMS, T=None, P=None):
    #main function defining basic parameterised TP-profile from 
    #PARAMS and INDEX. INDEX = [Chi, T,P]
    
        INDEX = self.TPindex
        
    
        if P is None:
            P = self.P
            
        if T is not None: 
            return T, P 

        if INDEX[1] > 1:
            P_params =  [self.params.tp_max_pres] + P_params + [np.min(P)]
            T_params = T_params + [T_params[-1]]
            #creating linear T-P profile
            T = np.interp(np.log(P[::-1]), np.log(P_params[::-1]), T_params[::-1])
            return T[::-1], P
        
        if INDEX[1] == 1:
            T = np.zeros_like(P)
            T += PARAMS[INDEX[0]+INDEX[1]]
            return T, P
        
        
        
        
        
        
        
        
        
            
    
    
    
    
    
        
