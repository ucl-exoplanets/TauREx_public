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
import library.library_emission as em
from library.library_general import *
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
        
        print 'z ',np.max(self.z)
        print 'p ',np.max(self.p)

        self.dzarray       = self.get_dz()

        if usedatagrid:
        #use wavelengthgrid of data or internal specgrid defined in data class
            self.set_lambdagrid(data['wavegrid'])
        else:
            self.set_lambdagrid(data['specgrid'])

    
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
        
        
    def path_integral(self, X=None, rho=None,temperature=None):
        
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.T#self.planet_temp

            
#         Tstar = 6000
        BB_star = self.F_star

        #constants 
        molnum = len(X[:,0])
        I_total    = np.zeros((self.nlambda))  
        tau        = np.zeros((self.nlayers,self.nlambda)) 
        dtau       = np.zeros((self.nlambda,self.nlambda)) 
        tau_total  = np.zeros((self.nlambda)) 
        
        #surface layer      
        BB_surf = em.black_body(self.specgrid,temperature[0])  
        sigma_array = self.get_sigma_array(temperature[0])
#            

        for k in range(self.nlayers):
                sigma_array = self.get_sigma_array(temperature[k])
                for i in range(molnum):
                    tau[0,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])
#                 print tau[0,0], X[i,0], rho[0], self.dzarray[0]
  
        exptau = np.exp(-1.0*tau[0,:])
        I_total += BB_surf*(exptau)
#         
#         print '-------------'
#         pl.figure(100)
#         pl.plot(self.specgrid,I_total)
#         pl.plot(self.specgrid,BB_surf,'r')
        

        
#         for j in xrange(self.nlayers-1,0,-1):
        for j in xrange(1,self.nlayers):

            for k in range(j,self.nlayers):
                sigma_array = self.get_sigma_array(temperature[k])
                for i in range(molnum):
#                     print rho[k]
                    tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])

#             sigma_array = self.get_sigma_array(temperature[j])
#             for i in range(molnum):
# #                     print rho[k]
#                     tau[j,:] += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j] )
                    
                
#                 pl.ion()
#                 pl.figure(100)
#                 pl.plot(rho[k],'x')
#                 pl.draw()
             
             
#         for j in range(self.nlayers-1):
            
#             dtau = tau[j+1,:] - tau[j,:]

            for i in range(molnum):
                dtau[j,:] += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j])

            
            exptau =  np.exp(-1.0*tau[j,:]) 
            
#             taulist.append(dtau)
#             print exptau
                
#             pl.ion()
#             pl.figure(100)
#             pl.plot(self.specgrid,exptau)
#             pl.draw()

             
            BB_layer = em.black_body(self.specgrid,temperature[j])           
            I_total += BB_layer*(exptau) * (dtau[j,:])
            
            
#             pl.ion()
#             pl.clf()
#             pl.figure(101)
#             pl.plot(self.specgrid,BB_layer*(exptau) * dtau[j,:],'b')
#             pl.plot(self.specgrid,I_total,'g')
#             pl.plot(self.specgrid,em.black_body(self.specgrid,temperature[0]),'k--')
#             pl.plot(self.specgrid,em.black_body(self.specgrid,temperature[-1]),'k-')
#             pl.plot(self.specgrid,BB_layer,'r')
#             pl.title(str(j))
#             pl.xscale('log')
#             pl.xlim([0.0,15.0])
#                   
# #             pl.figure(102)
# #             pl.plot(em.black_body(self.specgrid,temperature[0]),'k--')
# #             pl.plot(em.black_body(self.specgrid,temperature[-1]),'k-')
# #             pl.show()
#             pl.draw()
#                   
#             if j < 20:
#                 time.sleep(2)


        FpFs = (I_total/ BB_star) *(self.Rp/self.Rs)**2
        

        
        
#         ble_surf = (em.black_body(self.specgrid,temperature[0]) /BB_star) *(self.Rp/self.Rs)**2
#         ble_top = (em.black_body(self.specgrid,temperature[-1])/BB_star) *(self.Rp/self.Rs)**2
        
#         pl.ioff()
#         pl.figure(102)
#         pl.plot(self.specgrid,FpFs,label='spectrum')
#         pl.plot(self.specgrid,ble_surf,'k',label='bottom layer temp')
#         pl.plot(self.specgrid,ble_top,'k--',label='top layer temp')
#         pl.legend()
# #         pl.plot(self.specgrid,em.black_body(self.specgrid,1000))
# #         pl.plot(self.specgrid,sta,'r')
# #         pl.xscale('log')
# #      
# #         pl.figure(103)
# #         pl.plot(dtau[:,0])
# #         pl.title('dtau')
#         pl.show()

        return FpFs