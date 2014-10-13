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
from base import base    
import numpy as np
import pylab as pl
import library_emission as em
import library_general
from library_general import *
import time


class emission(base):

#initialisation
    def __init__(self,params,data,profile,usedatagrid=False):
        
        self.__ID__        = 'emission' #internal class identifier
        
        #loading data
        self.params        = params
        self.Rp            = params['planet_radius']
        self.Rs            = params['star_radius']

        self.n_gas         = data['ngas']
        self.specgrid      = data['specgrid']
        self.sigma_dict    = data['sigma_dict']
        self.X             = data['X']
        self.atmosphere    = data['atmosphere']
        self.F_star        = data['F_star']

        self.nlayers       = profile['nlayers']
        self.z             = profile['Z']        
        self.T             = profile['T']
        self.rho           = profile['rho']
        self.p             = profile['P']
        self.p_bar         = self.p * 1.0e-5 #convert pressure from Pa to bar
        
#         print 'z ',np.max(self.z)
#         print 'p ',np.max(self.p)

        self.dzarray       = self.get_dz()

        if usedatagrid:
        #use wavelengthgrid of data or internal specgrid defined in data class
            self.set_lambdagrid(data['wavegrid'])
        else:
            self.set_lambdagrid(data['specgrid'])
            
            
        #setting up static arrays for path_integral
        self.I_total    = np.zeros((self.nlambda))  
        self.tau        = np.zeros((self.nlayers,self.nlambda)) 
        self.dtau       = np.zeros((self.nlambda,self.nlambda)) 
        self.tau_total  = np.zeros((self.nlambda)) 

    
    #class methods 
        
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
        
#     @profile #line-by-line profiling decorator
    def path_integral(self, X=None, rho=None,temperature=None):
        
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.T#self.planet_temp

            

        BB_star = self.F_star

        #constants 
        I_total    = np.zeros((self.nlambda))  
        tau        = np.zeros((self.nlayers,self.nlambda)) 
        dtau       = np.zeros((self.nlambda,self.nlambda)) 
        tau_total  = np.zeros((self.nlambda)) 
        
        #surface layer      
        BB_surf = em.black_body(self.specgrid,temperature[0])  
        sigma_array = self.get_sigma_array(temperature[0])
#            

        for k in xrange(self.nlayers):
                sigma_array = self.get_sigma_array(temperature[k])
                for i in xrange(self.n_gas):
                    tau[0,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])

  
        exptau = np.exp(-1.0*tau[0,:])
        I_total += BB_surf*(exptau)


        for j in xrange(1,self.nlayers):

            for k in xrange(j,self.nlayers):
                sigma_array = self.get_sigma_array(temperature[k])
                for i in xrange(self.n_gas):
                    tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])


            for i in xrange(self.n_gas):
                dtau[j,:] += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j])

            
            exptau =  np.exp(-1.0*tau[j,:]) 
            
             
            BB_layer = em.black_body(self.specgrid,temperature[j])           
            I_total += BB_layer*(exptau) * (dtau[j,:])
            

        FpFs = (I_total/ BB_star) *(self.Rp/self.Rs)**2
 
        return FpFs