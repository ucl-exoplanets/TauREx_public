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
import ctypes as C
import library_emission as em
import library_general
from library_general import *
import time


class emission(object):

#initialisation
    def __init__(self,params,data,profile,usedatagrid=False):
        
        self.__ID__        = 'emission' #internal class identifier
        
        #type declations for arrays
        self.DTYPE = np.float
        
        #loading data
        self.params        = params
        self.Rp            = params['planet_radius']
        self.Rs            = params['star_radius']

        self.n_gas         = data['ngas']
        self.specgrid      = data['specgrid'].astype(self.DTYPE)
        self.sigma_dict    = data['sigma_dict']
        self.X             = data['X'].astype(self.DTYPE)
        self.atmosphere    = data['atmosphere']
        self.F_star        = data['F_star'].astype(self.DTYPE)

        self.nlayers       = profile['nlayers']
        self.z             = profile['Z'].astype(self.DTYPE)        
        self.T             = profile['T'].astype(self.DTYPE)
        self.rho           = profile['rho'].astype(self.DTYPE)
        self.p             = profile['P'].astype(self.DTYPE)
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
        self.I_total    = np.zeros((self.nlambda),dtype=self.DTYPE)  
        self.tau        = np.zeros((self.nlayers,self.nlambda),dtype=self.DTYPE) 
        self.dtau       = np.zeros((self.nlambda,self.nlambda),dtype=self.DTYPE) 
        self.tau_total  = np.zeros((self.nlambda),dtype=self.DTYPE) 

        #loading c++ pathintegral library for faster computation
        if params.trans_cpp:
            # self.cpathlib = C.cdll.LoadLibrary('./library/pathintegral_test.so')
            self.cpathlib = C.CDLL('./library/pathintegral_emission.so',mode=C.RTLD_GLOBAL)
            self.sigma_array_c, self.sig_tempgrid = self.get_sigma_array_c()
    
    #class methods 
        
    def set_lambdagrid(self,GRID):
    #sets internal memory of wavelength grid to be used
        self.lambdagrid = GRID
        self.nlambda = len(GRID)       
        
    
    def get_sigma_array(self,temperature):
    #getting sigma array from sigma_dic for given temperature 
#         print temperature 
        return self.sigma_dict[find_nearest(self.sigma_dict['tempgrid'],temperature)[0]]
    
    def get_sigma_array_c(self):
    #generating 3D sigma_array from sigma_dict for c++ path integral
        tempgrid = self.sigma_dict['tempgrid']
        OUT = np.zeros((len(tempgrid),2,len(self.specgrid)),dtype=np.float64)
        
        c=0
        for t in tempgrid:
            OUT[c,:,:] = self.sigma_dict[t]
            c += 1
        
        return OUT, np.asarray(tempgrid,dtype=np.float64)
    
        
    def get_dz(self):
        
        dz = []
        for i in range(len(self.z)-1):
            dz.append(self.z[i+1]-self.z[i])
            
        dz.append(dz[-1])
        
        return np.asarray(dz,dtype=self.DTYPE)
        
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
        self.I_total[:]   = 0
        self.tau[:]       = 0
        self.dtau[:]      = 0
        self.tau_total[:] = 0
        
        #surface layer      
        BB_surf = em.black_body(self.specgrid,temperature[0])  
        sigma_array = self.get_sigma_array(temperature[0])

        for k in xrange(self.nlayers):
                if temperature[k] != temperature[k-1]:
                    sigma_array = self.get_sigma_array(temperature[k])
                    
                for i in xrange(self.n_gas):
                    self.tau[0,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])
                    
        exptau = np.exp(-1.0*self.tau[0,:])
        self.I_total += BB_surf*(exptau)

#         for i in range(100):
#             print BB_surf[i], self.tau[0,i], np.exp(-1.0*self.tau[0,i])
            
        #other layers
        BB_layer = BB_surf
        sigma_array = self.get_sigma_array(temperature[0])
        for j in xrange(1,self.nlayers):

            for k in xrange(j,self.nlayers):
                if temperature[k] != temperature[k-1]:
                    sigma_array = self.get_sigma_array(temperature[k])
                    
                for i in xrange(self.n_gas):
                    self.tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])
                if k is j:
                    self.dtau[j,:] = self.tau[j,:]
                      
#             for i in xrange(self.n_gas):
#                 self.dtau[j,:] += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j])

            
            exptau =  np.exp(-1.0*self.tau[j,:]) 
            
            if temperature[j] != temperature[j-1]: 
                BB_layer = em.black_body(self.specgrid,temperature[j])           
            self.I_total += BB_layer*(exptau) * (self.dtau[j,:])
            

        FpFs = (self.I_total/ BB_star) *(self.Rp/self.Rs)**2
 
        return FpFs
    
    def cpath_integral(self, X = None, rho = None, temperature= None):
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.T

        
        #casting changing arrays to c++ pointers
        Xs1,Xs2 = shape(X)
        Xnew = zeros((Xs1+1,Xs2))
        Xnew[:-1,:] = X
        cX = cast2cpp(Xnew)
        crho = cast2cpp(rho)
        ctemperature = cast2cpp(temperature)
        cF_star =cast2cpp(self.F_star)
        cspecgrid = cast2cpp(self.specgrid)
        #casting fixed arrays and variables to c++ pointers
        
#         csigma_array = cast2cpp(self.sigma_array_c)
        csigma_array = self.sigma_array_c.ctypes.data_as(C.POINTER(C.c_double))

        csig_tempgrid = cast2cpp(self.sig_tempgrid)
        cdzarray = cast2cpp(float64(self.dzarray))
        znew = zeros((len(self.z)))
        znew[:] = self.z
        cz  = cast2cpp(znew)
        cdz  = cast2cpp(self.dzarray)
#         cRsig = cast2cpp(self.Rsig)
#         cCsig = cast2cpp(self.Csig)
            
        cRp = C.c_double(self.Rp)
        cRs = C.c_double(self.Rs)
        clinecount = C.c_int(self.nlambda)
        cnlayers = C.c_int(self.nlayers)
        cn_gas = C.c_int(len(X[:,0]))
        cn_lambda = C.c_int(len(self.specgrid))
        cn_sig_temp= C.c_int(len(self.sig_tempgrid))
        
        #setting up output array
        FpFs = zeros((self.nlambda),dtype=float64)

        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int_emission
        
        cpath_int(cX,crho,ctemperature,cF_star,cspecgrid,csigma_array,cdzarray,
                  cn_lambda,cRp,cRs,cnlayers,cn_gas,csig_tempgrid,cn_sig_temp, C.c_void_p(FpFs.ctypes.data))
        
        OUT = zeros((self.nlambda))
        OUT[:] = FpFs
        del(FpFs)

        return OUT
        
