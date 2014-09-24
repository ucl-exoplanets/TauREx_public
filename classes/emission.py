################################################
#class emission
#
# Calculates the emission radiative transfer 
# forward model
#
# Input: -
#
#
# Output: -  
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013 
#
################################################

#loading libraries     
import numpy as np
import pylab as pl
import library.library_emission as em
from library.library_general import *


class emission(object):

#initialisation
    def __init__(self,params,data,profile,usedatagrid=False):
        #loading data
        self.params        = params
        self.Rp            = params['planet_radius']
        self.Rs            = params['star_radius']

        self.n_gas         = data['ngas']
        self.specgrid      = data['specgrid']
        self.sigma_dict    = data['sigma_dict']
        self.X             = data['X']
        self.atmosphere    = data['atmosphere']

        self.nlayers       = profile['nlayers']
        self.z             = profile['Z']        
        self.t             = profile['T']
        self.rho           = profile['rho']
        self.p             = profile['P']
        self.p_bar         = self.p * 1.0e-5 #convert pressure from Pa to bar
        
        print 'z ',np.max(self.z)
        print 'p ',np.max(self.p)

        self.dzarray       = self.get_dz()

        if usedatagrid:
        #use wavelengthgrid of data or internal specgrid defined in data class
            self.set_lambdagrid(data['wavegrid'])
        else:
            self.set_lambdagrid(data['specgrid'])

#basic class methods and overloading
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

    def reset(self,data):
    #allows to reset the original instance to reflect changes in the data instance
    #this avoids an initialisation of a separate instance.
    # noinspection PyArgumentList
        self.__init__(self.params,data)
    
    
    def set_lambdagrid(self,GRID):
    #sets internal memory of wavelength grid to be used
        self.lambdagrid = GRID
        self.nlambda = len(GRID)       
        
    
    def get_sigma_array(self,temperature):
    #getting sigma array from sigma_dic for given temperature 
#         print temperature 
        return self.sigma_dict[find_nearest(self.sigma_dict['tempgrid'],temperature)[0]]
        
        
    def get_dz(self):
        
        dz = []
        for i in range(len(self.z)-1):
            dz.append(self.z[i+1]-self.z[i])
            
        dz.append(dz[-1])
        
        return np.asarray(dz)
        
        
    def path_integral(self, X=None, rho=None,temperature=None):
        
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.t#self.planet_temp
            
        pl.figure(100)
        pl.plot(self.z,rho)
        
        pl.figure(101)
        pl.plot(self.z,self.t)

            
        Tstar = 4500
        BB_star = em.black_body(self.specgrid,Tstar)

        #constants 
        molnum = len(X[:,0])
        I_total    = np.zeros((self.nlambda))  
        tau        = np.zeros((self.nlambda)) 
        tau_total  = np.zeros((self.nlambda)) 
        

        for j in range(self.nlayers):
            tau[:] = 0.0
            sigma_array = self.get_sigma_array(temperature[j]) #selecting correct sigma_array for temperature
            for i in range(molnum):
                tau += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j])

            exptau =  np.exp(-1.0*tau)    
            BB_layer = em.black_body(self.specgrid,temperature[j])
            
            I_total += BB_layer*exptau
        
        
        FpFs = (self.Rp**2 *I_total) / (self.Rs**2 * BB_star)

        
#         pl.figure(101)
#         pl.plot(self.specgrid,self.Rp**2 *I_total)
#         pl.plot(self.specgrid,self.Rs**2 * BB_star,'r')
        
        pl.figure(102)
        pl.plot(self.specgrid,FpFs)
        pl.show()
        
        
        return 1