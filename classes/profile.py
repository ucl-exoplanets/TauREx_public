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
        self.pta    = data.pta
        self.X      = data.X
        
        if params.tp_var_atm == False:
            self.nlayers = data.nlayers
            self.pta     = data.pta
            self.X       = data.X
        
        
                
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
        
