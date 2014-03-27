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
import numpy
from numpy import *



class profile(object):


    def __init__(self,params,data):
        
        self.params = params
        
        #constants
        self.botlzman = 1.3806488e-23# m2 kg s-2 K-1
        
        #derived values
        self.scaleheight = self.get_scaleheight(params.planet_temp, params.planet_grav, params.planet_mu)
     
        
        
        if params.tp_var_atm == True:
            self.nlayers = int(params.tp_atm_levels)
            self.ngas    = int(params.tp_num_gas)
            self.pta     = self.setup_pta_grid()
            self.P       = self.pta[:,0]
            self.T       = self.pta[:,1]
            self.Z       = self.pta[:,2]
            self.X       = zeros((self.nlayers,self.ngas))
            self.X      += 1e-5  #setting up initial mixing ratios
            
        else:           
            self.nlayers = data.nlayers
            self.pta     = data.pta
            self.T       = data.pta[:,1]
            self.X       = data.X       
            self.P       = self.pta[:,0]
        
                
#basic class methods and overloading
    def list(self,name=None):
        if name==None:
            return dir(self)[2:-1]
        else:
            lst = dir(self)
            return filter(lambda k: name in k, lst)
        
    def __getattribute__(self,name):
        return object.__getattribute__(self, name)
    
    def __getitem__(self,name):
        return self.__dict__[name] 

    def reset(self,data):
    #allows to reset the original instance to reflect changes in the data instance
    #this avoids an initialisation of a separate instance.
        self.__init__(self.params,data)
        
    def get_scaleheight(self,T_aver,surf_g,mmw):
        
        return (self.botlzman*T_aver)/(mmw*surf_g)
        
        
        
    def setup_pta_grid(self):
        
        MAX_P    = self.params.tp_max_pres
        N_SCALE  = self.params.tp_num_scale
        N_LAYERS = self.nlayers
        
        max_z    = N_SCALE * self.scaleheight
#         dz       = max_z / N_LAYERS
        
        #generatinng altitude-pressure array
        PTA_arr = zeros((N_LAYERS,3))
        PTA_arr[:,2] = linspace(0,max_z,num=N_LAYERS)
        PTA_arr[:,0] = MAX_P * exp(-PTA_arr[:,2]/self.scaleheight)
        PTA_arr[:,1] = self.params.planet_temp

        
        return PTA_arr
        
        
                
        
    def get_rho(self,T=None,P=None):
        
        if P == None:
            P = self.P
        if T == None:
            T = self.params.planet_temp
            
        return  (P)/(self.botlzman*T)   
        
        
        
