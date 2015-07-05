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
#loading classes
from base import base    
import numpy as np
import pylab as pl
import ctypes as C
import library_emission as em
import library_general as gen

try:
    import cy_pathintegral_emission
    from cy_pathintegral_emission import path_integral as cython_path_integral
    cythonised = True
except ImportError:
    cythonised = False

# from library_general import *
import time
import logging

class emission(base):

#initialisation
#    def __init__(self,params,data,profile,usedatagrid=False):
    def __init__(self, atmosphere, data=None, params=None, usedatagrid=False):

        logging.info('Initialise object emission')

        if params:
            self.params = params
        else:
            self.params = atmosphere.params #get params object from profile

        if data:
            self.data = data
        else:
            self.data = atmosphere.data    # get data object from profile

        self.__ID__        = 'emission' #internal class identifier

        #type declations for arrays
        self.DTYPE = np.float

        #loading data
        self.atmosphere    = atmosphere
        self.Rp            = self.params.planet_radius
        self.Rs            = self.params.star_radius

        self.n_gas         = self.atmosphere.ngas
        self.specgrid      = self.data.specgrid.astype(self.DTYPE)
        self.sigma_dict    = self.data.sigma_dict
        self.X             = self.atmosphere.X.astype(self.DTYPE)
        self.F_star        = self.data.F_star.astype(self.DTYPE)
        
        self.nlayers       = self.atmosphere.nlayers
        self.z             = self.atmosphere.z.astype(self.DTYPE)
        self.T             = self.atmosphere.T.astype(self.DTYPE)
        self.rho           = self.atmosphere.rho.astype(self.DTYPE)
        self.P             = self.atmosphere.P.astype(self.DTYPE)
        self.P_bar         = self.P * 1.0e-5 #convert pressure from Pa to bar

        self.dzarray       = self.get_dz()

        if usedatagrid:
        #use wavelengthgrid of data or internal specgrid defined in data class
            self.set_lambdagrid(self.data.wavegrid)
        else:
            self.set_lambdagrid(self.data.specgrid)


        #setting up static arrays for path_integral
        self.I_total    = np.zeros((self.nlambda),dtype=self.DTYPE)
        self.tau        = np.zeros((self.nlayers,self.nlambda),dtype=self.DTYPE)
        self.dtau       = np.zeros((self.nlayers,self.nlambda),dtype=self.DTYPE)
        self.tau_total  = np.zeros((self.nlayers,self.nlambda),dtype=self.DTYPE)

        #loading c++ pathintegral library for faster computation
        if self.params.trans_cpp:
            # self.cpathlib = C.cdll.LoadLibrary('./library/pathintegral_test.so')
            self.cpathlib = C.CDLL('./library/pathintegral_emission.so',mode=C.RTLD_GLOBAL)
            self.sigma_array_c, self.sig_tempgrid = self.get_sigma_array_c()

        # set forward model function
        if cythonised:
            self.model = self.cy_path_integral
        else:
            self.model = self.path_integral


    #class methods
    def set_lambdagrid(self,GRID):
        #sets internal memory of wavelength grid to be used
        self.lambdagrid = GRID
        self.nlambda = len(GRID)


    def get_sigma_array(self,temperature):
    #getting sigma array from sigma_dic for given temperature
        return self.sigma_dict[gen.find_nearest(self.sigma_dict['tempgrid'],temperature)[0]]

    def get_sigma_array_c(self):
    #generating 3D sigma_array from sigma_dict for c++ path integral
        tempgrid = self.sigma_dict['tempgrid']
        OUT = np.zeros((len(tempgrid),self.n_gas,len(self.specgrid)),dtype=np.float64)
        
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
    
    def get_contribution_function(self):
        '''
        Function re-running forward model with current setup and
        explicitely returning the path-integral byproducts: tau, tau_total, dtau
        '''    
        model = self.path_integral() #output not needed but needs running to set arrays
    
        return self.tau, self.tau_total, self.dtau
        
        
        
#     @profile #line-by-line profiling decorator
    def path_integral(self, X=None, rho=None,temperature=None):
        
        if X is None:
            X = self.atmosphere.X
        if rho is None:
            rho = self.atmosphere.rho
        if temperature is None:
            temperature = self.atmosphere.T#self.planet_temp


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
            
        #other layers
        BB_layer = BB_surf
        
        sigma_array = self.get_sigma_array(temperature[0])
        for j in xrange(1,self.nlayers):   
            if temperature[j] != temperature[j-1]: 
                BB_layer = em.black_body(self.specgrid,temperature[j]) 
                            
            for k in xrange(j,self.nlayers):
                if temperature[k] != temperature[k-1]:
                    sigma_array = self.get_sigma_array(temperature[k])  
                    
                if j is k:
                    for i in xrange(self.n_gas):
                        self.tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])
                        self.dtau[j,:] += (sigma_array[i,:] * X[i,j] * rho[j] * self.dzarray[j])
                else:
                    for i in xrange(self.n_gas):
                        self.tau[j,:] += (sigma_array[i,:] * X[i,k] * rho[k] * self.dzarray[k])  
            
            exptau =  np.exp(-1.0*self.tau[j,:]) 
             
              
            
            self.tau_total[j,:] = BB_layer*(exptau) * (self.dtau[j,:])
            self.I_total += BB_layer*(exptau) * (self.dtau[j,:])


        self.FpFs = (self.I_total/ BB_star) *(self.Rp/self.Rs)**2
        return self.FpFs
    
    
    
    def cy_path_integral(self, X=None, rho=None,temperature=None):
        
        if X is None:
            X = self.atmosphere.X
        if rho is None:
            rho = self.atmosphere.rho
        if temperature is None:
            temperature = self.atmosphere.T#self.planet_temp
            
        I_total, tau_total = cython_path_integral(X, rho, temperature, self.sigma_dict, self.F_star, self.specgrid, self.dzarray,self.nlambda,self.nlayers, self.n_gas)
        
        self.FpFs = (I_total/ self.F_star) *(self.Rp/self.Rs)**2
        return self.FpFs
    
    
    def cpath_integral(self, X = None, rho = None, temperature= None):
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.T

        
        #casting changing arrays to c++ pointers
        Xs1,Xs2 = np.shape(X)
        Xnew = np.zeros((Xs1+1,Xs2))
        Xnew[:-1,:] = X
        cX = gen.cast2cpp(Xnew)
        crho = gen.cast2cpp(rho)
        ctemperature = gen.cast2cpp(temperature)
        cF_star = gen.cast2cpp(self.F_star)
        cspecgrid = gen.cast2cpp(self.specgrid)
        #casting fixed arrays and variables to c++ pointers
        
#         csigma_array = cast2cpp(self.sigma_array_c)
        csigma_array = self.sigma_array_c.ctypes.data_as(C.POINTER(C.c_double))

        csig_tempgrid = gen.cast2cpp(self.sig_tempgrid)
        cdzarray = gen.cast2cpp(np.float64(self.dzarray))
        znew = np.zeros((len(self.z)))
        znew[:] = self.z
        cz  = gen.cast2cpp(znew)
        cdz  = gen.cast2cpp(self.dzarray)
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
        FpFs = np.zeros((self.nlambda),dtype=np.float64)

        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int_emission
        
        cpath_int(cX,crho,ctemperature,cF_star,cspecgrid,csigma_array,cdzarray,
                  cn_lambda,cRp,cRs,cnlayers,cn_gas,csig_tempgrid,cn_sig_temp, C.c_void_p(FpFs.ctypes.data))
        
        OUT = np.zeros((self.nlambda))
        OUT[:] = FpFs
        del(FpFs)

        return OUT
        
