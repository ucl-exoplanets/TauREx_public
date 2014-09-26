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
            temperature = self.T#self.planet_temp
            
#         pl.figure(100)
#         pl.plot(self.z,rho)
#         
#         pl.figure(101)
#         pl.plot(self.z,self.t)

            
        Tstar = 6000
        BB_star = em.black_body(self.specgrid,Tstar)

        #constants 
        molnum = len(X[:,0])
        I_total    = np.zeros((self.nlambda))  
        tau        = np.zeros((self.nlayers+1,self.nlambda)) 
        tau_total  = np.zeros((self.nlambda)) 
        
        #surface layer
        
        BB_surf = em.black_body(self.specgrid,temperature[0])  
        sigma_array = sigma_array = self.get_sigma_array(temperature[0])
           
        for i in range(molnum):
                tau[0,:] += (sigma_array[i,:] * X[i,0] * rho[0] * self.dzarray[0])
        exptau = np.exp(-1.0*tau[0,:])
        I_total += BB_surf*exptau
        
        for j in range(1,self.nlayers):
#             tau[:] = 0.0
            

            for k in range(j,self.nlayers):
                sigma_array = self.get_sigma_array(temperature[k]) #selecting correct sigma_array for temperature
                for i in range(molnum):
#                     print rho[k]
                    tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k] )
                    
                
#                 pl.ion()
#                 pl.figure(100)
#                 pl.plot(rho[k],'x')
#                 pl.draw()
             
             
        for j in range(1,self.nlayers):
            exptau = 0.0
            
            dtau = tau[j,:] - tau[j-1,:]
            print dtau
            
            exptau =  np.exp(-1.0*tau[j,:]) 
            
#             taulist.append(dtau)
#             print exptau
                
#             pl.ion()
#             pl.figure(100)
#             pl.plot(self.specgrid,exptau)
#             pl.draw()

             
            BB_layer = em.black_body(self.specgrid,temperature[j])
             
            I_total += BB_layer*(-1.0*exptau) * dtau
            

        
#         pl.show()
        FpFs = (I_total/ BB_star) *(self.Rp/self.Rs)**2
        
        
#         pl.figure(101)
#         pl.plot(self.specgrid,self.Rp**2 *I_total)
#         pl.plot(self.specgrid,self.Rs**2 * BB_star,'r')
        
        
        pla = em.black_body(self.specgrid,1400) 
        sta = em.black_body(self.specgrid,6000) 
        
        pla2 = em.black_body(self.specgrid,3000)
        
        ble = pla/sta *(self.Rp/self.Rs)**2
        ble2 = pla2/sta *(self.Rp/self.Rs)**2
        
        ble_surf = BB_surf/sta *(self.Rp/self.Rs)**2
        
        
        pl.figure(102)
        pl.plot(self.specgrid,FpFs)
        pl.plot(self.specgrid,ble,'k')
        pl.plot(self.specgrid,ble2,'k--')
#         pl.plot(self.specgrid,em.black_body(self.specgrid,1000))
#         pl.plot(self.specgrid,sta,'r')
#         pl.xscale('log')
     
#         pl.figure(103)
#         pl.plot(taulist,'x')
#         pl.title('tau')
        pl.show()
        
        
        
        return 1