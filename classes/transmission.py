################################################
#class transmission
#
# Calculates the transmission radiative transfer 
# forward model. Based on tau.cpp (Morgan Hollis).
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
import numpy
import copy
from numpy import *
from pylab import *
import itertools
import time
import ctypes as C
import library_general,library_transmission
from library_general import *
from library_transmission import *
import logging

# from mpi4py import MPI

class transmission(base):

#    def __init__(self,params,data,profile,usedatagrid=False):
    def __init__(self, profile, data=None, params=None, usedatagrid=False):

        logging.info('Initialise object transmission')

        self.__ID__        = 'transmission' #internal class identifier

        if params:
            self.params = params
        else:
            self.params = profile.params #get params object from profile

        if data:
            self.data = data
        else:
            self.data = profile.data    # get data object from profile

        self.Rp            = self.params.planet_radius   # @todo As long as they don't change, better not to reassign params vars to new local variables
        self.Rs            = self.params.star_radius     #

        self.ngas          = self.data.ngas               # @todo If var comes from data object, better to be explicit later as well
        self.sigma_dict    = self.data.sigma_dict
        self.X             = self.data.X
        self.atmosphere    = self.data.atmosphere
        
        self.nlayers       = profile.nlayers              # @todo same here
        self.z             = profile.z
        self.rho           = profile.rho
        self.P             = profile.P
        self.P_bar         = self.P * 1.0e-5 #convert pressure from Pa to bar # @todo move to tp_profile
        
        self.planet_temp   = int(self.params.planet_temp)

        #cloud specific parameters.
        self.include_cld   = self.params.in_include_cld # @todo and here
        if self.include_cld:
            self.cld_m         = self.params.in_cld_m #refractive index for cloud cross sections
            self.cld_a         = self.params.in_cld_a #cloud particle size (microns)
            self.cld_upbound   = self.params.in_cld_pressure[1] #pressure (in bar) of upper pressure/lower altitude cloud bound
            self.cld_lowbound  = self.params.in_cld_pressure[0] #pressure (in bar) of lower pressure/upper altitude cloud bound
        
        if usedatagrid:
            #use wavelengthgrid of data or internal specgrid defined in data class
            self.set_lambdagrid(self.data.wavegrid)
        else:
            self.set_lambdagrid(self.data.specgrid)

        #calculating optical path lengths
        self.dlarray,self.iteridx = self.get_path_length()
        self.dz            = self.get_dz()

        #calculating rayleigh scattering cross sections
        self.mol_list      = self.atmosphere['mol'].keys()
        self.ray_thres     = 50.0 * (self.atmosphere['mol'][self.mol_list[0]]['radius'] * 1.0e6);
        self.Rsig          = self.get_Rsig()

        #calculating collision induced absorption cross sections
        if self.params.in_include_cia:
            self.cia       = self.data.cia
            self.Csig      = self.get_Csig()
        else:
            self.Csig      = zeros((self.nlambda))
            
        #calculating cloud cross sections
        if self.include_cld:
            self.Cld_sig = self.get_cloud_sig()
        else:
           self.Cld_sig = zeros((self.nlambda)) #unused but needed to cast to cpp code
        
        #loading c++ pathintegral library for faster computation
        if self.params.trans_cpp:
            # self.cpathlib = C.cdll.LoadLibrary('./library/pathintegral_test.so')
            self.cpathlib = C.CDLL('./library/pathintegral.so', mode=C.RTLD_GLOBAL)


#class methods
    def set_lambdagrid(self, grid):
    #sets internal memory of wavelength grid to be used
        self.lambdagrid = grid
        self.nlambda = len(grid)

    def get_path_length(self):
        #calculates the layerlength 
        dlarray = []
        jl = []
        kl = []
        cl = []
        count = 0
        for j in range(self.nlayers):    # loop through atmosphere layers
            for k in range(1,self.nlayers-j):    # loop through each layer to sum up path length
                dl = 2.0 * (sqrt(pow((self.Rp + self.z[k+j]),2) - pow((self.Rp + self.z[j]),2)) -
                            sqrt(pow((self.Rp + self.z[k-1+j]),2) - pow((self.Rp + self.z[j]),2)))
                dlarray.append(dl)
                jl.append(j)
                kl.append(k)
                cl.append(count)
                count += 1
         
        iteridx = zip(cl,jl,kl)
        dlarray = asarray(dlarray)      
#         dlarray[isnan(dlarray)] = 0


        return dlarray, iteridx
    
    def get_dz(self):
        #calculate area of circles of atmos
        dz = zeros((self.nlayers))
        for j in range(self.nlayers-1):
            dz[j] = self.z[j+1] - self.z[j]
        return dz
    
    def get_Rsig(self):
    #calculating rayleigh scattering cross-sections
        Rsig = zeros((self.nlambda))
        count = 0
        for wl in self.lambdagrid:
            if wl < self.ray_thres: 
                pass
            else:
                for mol in self.mol_list:
                    Rsig[count] += ( self.atmosphere['mol'][mol]['frac'] * self.scatterRayleigh(wl,mol) )
#             print wl,Rsig[count]
            count += 1
        
        return Rsig

    def get_Csig(self):
    #calculating collisional induced absorption
        return self.scatterCIA(self.cia[:,1],self.atmosphere['mol']['H2']['frac'])      
    
    def get_cloud_sig(self):
    #calculating could cross sections
        return (128.0* np.pi**5 * self.cld_a**6) / (3.0 * self.lambdagrid**4) * ((self.cld_m**2 -1.0)/(self.cld_m**2 +2.0))**2

        

    def get_sigma_array(self,temperature):
    #getting sigma array from sigma_dic for given temperature 
#         print temperature 
        return self.sigma_dict[find_nearest(self.sigma_dict['tempgrid'],temperature)[0]]
        

    def path_integral(self, X = None, rho = None, temperature=None):

        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.planet_temp
            
        #selecting correct sigma_array for temperature
        sigma_array = self.get_sigma_array(temperature)
        
        #setting up arrays
        absorption = zeros((self.nlambda))
        tau        = zeros((self.nlayers))  
        Rtau       = zeros((self.nlayers))
        Ctau       = zeros((self.nlayers))
        cld_tau    = zeros((self.nlayers))
        exptau     = zeros((self.nlayers))

        molnum = len(X[:,0])
        #running loop over wavelengths     
        for wl in range(self.nlambda):
            tau[:]     = 0.0 
            exptau[:]  = 0.0
            Rtau[:]    = 0.0
            cld_tau[:] = 0.0
         
            for c,j,k in self.iteridx:
                #optical depth due to gasses
                for i in range(molnum):
                    tau[j] += (sigma_array[i,wl] * X[i,j+k] * rho[j+k] * self.dlarray[c])
                
                Rtau[j] = self.Rsig[wl] * rho[j+k] * self.dlarray[c] #optical depth due to Rayleigh scattering
                Ctau[j] = self.Csig[wl] * (rho[j+k] **2) * self.dlarray[c]  # calculating CIA optical depth
                
                #Calculating cloud opacities, single cloud layer implementation only. can be upgraded or ... not
                if self.include_cld:                    
                    if self.P_bar[k+j] < self.cld_upbound and self.P_bar[k+j] > self.cld_lowbound:
                        #= log(cloud density), assuming linear decrease with decreasing log pressure
                        #following Ackerman & Marley (2001), Fig. 6
                        cld_log_rho = interp_value(np.log(self.P_bar[k+j]),np.log(self.cld_lowbound),np.log(self.cld_upbound),-6.0,-1.0) 
                        cld_tau[j] += ( self.Cld_sig[wl] * (self.dlarray[c]*1.0e2) * (exp(cld_log_rho)*1.0e-6) )    # convert path lenth from m to cm, and density from g m^-3 to g cm^-3
                
                #adding all taus together
                tau[j] += Rtau[j]
                tau[j] += Ctau[j]
                tau[j] += cld_tau[j]
 
                exptau[j]= exp(-tau[j])
  
            integral = 2.0* sum(((self.Rp+self.z)*(1.0-exptau)*self.dz))
#             print self.z, self.dz, exptau
            
            absorption[wl] = (self.Rp**2 + integral) / (self.Rs**2)

        return absorption
        
    # @profile #line-by-line profiling decorator
    def cpath_integral(self, X = None, rho = None, temperature=None):
        if X is None:
            X = self.X
        if rho is None:
            rho = self.rho
        if temperature is None:
            temperature = self.planet_temp
            
        #selecting correct sigma_array for temperature
        sigma_array = self.get_sigma_array(temperature)
        
        #casting changing arrays to c++ pointers
        Xs1,Xs2 = shape(X)
        Xnew = zeros((Xs1+1,Xs2))
        Xnew[:-1,:] = X
        cX = cast2cpp(Xnew)
        crho = cast2cpp(rho)

        #casting fixed arrays and variables to c++ pointers
        csigma_array = cast2cpp(sigma_array)
        cdlarray = cast2cpp(self.dlarray)
        znew = zeros((len(self.z)))
        znew[:] = self.z
        cz   = cast2cpp(znew)
        cdz  = cast2cpp(self.dz)
        cRsig = cast2cpp(self.Rsig)
        cCsig = cast2cpp(self.Csig)
            
        cRp = C.c_double(self.Rp)
        cRs = C.c_double(self.Rs)
        clinecount = C.c_int(self.nlambda)
        cnlayers = C.c_int(self.nlayers)
        cn_gas = C.c_int(len(X[:,0]))

        #setting and casting cloud paramters
        cP_bar    = cast2cpp(self.P_bar)
        cCld_sig  = cast2cpp(self.Cld_sig)
        if self.include_cld:
            cInclude_cld  = C.c_int(1)
            cCld_lowbound = C.c_double(self.cld_lowbound)
            cCld_upbound  = C.c_double(self.cld_upbound)
        else:
            cInclude_cld  = C.c_int(0)
            cCld_lowbound = C.c_double(0)
            cCld_upbound  = C.c_double(0)

        #setting up output array
        absorption = zeros((self.nlambda),dtype=float64)

        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int

        #running c++ path integral
        cpath_int(csigma_array,cdlarray,cz,cdz,cRsig,cCsig,cX,crho,cRp,cRs,\
                  clinecount,cnlayers,cn_gas,cInclude_cld,cCld_lowbound,\
                  cCld_upbound,cP_bar,cCld_sig,C.c_void_p(absorption.ctypes.data))

        # del(cX); del(crho); del(csigma_array); del(cz); del(cRsig); del(cCsig); del(cRp); del(cRs);
        # del(clinecount); del(cnlayers); del(cn_gas); del(cpath_int)

        OUT = zeros((self.nlambda))
        OUT[:] = absorption
        del(absorption)

        return OUT
        

    def scatterRayleigh(self,wl,NAME):
        '''
         Formula from Liou 2002, 'An Introduction to Atmospheric Radiation', pp.92-93. Also uses 'minimum volume' approximation pg.97, 
            N_dens = 1 / V_particle .
        Optical depth given by tau = sigma * L * c ; sigma = abs cross-section (m^2), L = path length (m), c = concentration (m^-3)
        NB This is for bulk atmos scattering ONLY (assumptions: particles much smaller than wavelength, gas sufficiently dense), 
        cloud Rayleigh + Mie included in scatterMie function. 
        IN: Wavelength (in um), path length (in m)
        OUT: Rayleigh scattering opacity cross-section per particle (in m^2)
        '''
        #sigma_R=0.0                 # Rayleigh absorption coefficient (from Liou, An Introduction to Atmospheric Radiation)
        wl *= 1.0e-6  #convert wavelengths to m
        rad=self.atmosphere['mol'][NAME]['radius']    #molecular radius (m)
        r_ind=self.atmosphere['mol'][NAME]['ridx']    #molecular refractive index
        r_sq=r_ind**2 
        r_red = (r_sq-1) / (r_sq+2);

        delta = 0.035;               #molecular anisotropy factor
        f_delta = (6.0+(3.0*delta)) / (6.0-(7.0*delta));   #King correction factor
        if( NAME == 'He'): f_delta = 1.0;                  #no asymmetry for helium molecules

        # Find cross-section 
        sigma_R = (128.0/3.0) * (pow(pi,5) * pow(rad,6) / pow(wl,4)) *r_red*r_red *f_delta    #gives sigma_R in m^2
        
        return sigma_R
        
        
    def scatterCIA(self,coeff,amount):
        '''
        Optical depth given by tau = alpha * L * c_1 * c_2 ; alpha = abs coeff (cm^5 mol^-2), L = path length (cm), c_i = concentration of collider i (mol cm^-3)
        IN: CIA coeffs in (cm^-1 amagat^-2), grid wavelength (in um), path length (in m) and total number density dz (in m^-3)
        OUT: H2-H2 collision-induced absorption coefficient (in m^5 mol^-2)
        '''
        AMA=2.68676e+25 #Amagat (molecules m^-3)
        conv_factor = 1.0 / pow((AMA*1.0e-6),2) #converting from (cm^-1 amagat^-2) to (cm^5 mol^-2)...
        alpha = coeff * conv_factor
        alpha *= (amount**2) * 1.0e-10 #e.g. composition 85% H2, and convert from cm^5 to m^5
        return alpha
    
    
###################################
#old code sniplets do not delete

            
#             for j in range(self.nlayers-1):    
#                 for k in range(1,self.nlayers-j):
#                     # Calculate half-path length, and double (from system geometry) to get full path distance
# #                     dl = 2.0 * (sqrt(pow((self.Rp + self.z[k+j]),2) - pow((self.Rp + self.z[j]),2)) - sqrt(pow((self.Rp + self.z[k-1+j]),2) - pow((self.Rp + self.z[j]),2)))
#                      
# #                     print dl , self.dlarray[k+j]
# #                     for l in range(self.ngas):
#  
#  
#                     tau[j] += (self.sigma_array[:,wl] * transpose(self.X[k+j,:]) * self.rho[k+j] * self.dlarray[cidx])
#                     cidx += 1
#                 exptau[j]= exp(-tau[j])
        
            #calculate area of circles of atmos
#             dz = zeros((self.nlayers))
#             for j in range(self.nlayers-1):
#                 dz[j] = self.z[j+1] - self.z[j]
#             dz[-1] = dz[-2]

#             integral = 0.0
#             for j in range(self.nlayers):
#                 integral += ((self.Rp+self.z[j])*(1.0-exptau[j])*self.dz[j])
#             integral *= 2.0
            
            
            
#         for wl in range(self.nwave):     
#             tau[:] = 0.0 
#             exptau[:] = 0.0
# 
#             for c,j,k in self.iteridx:
#                 tau[j] += (self.sigma_array[:,wl] * transpose(self.X[j+k,:]) * self.rho[j+k] * self.dlarray[c])
#                 exptau[j]= exp(-tau[j])
# 
#             integral = 2.0* sum(((self.Rp+self.z)*(1.0-exptau)*self.dz))
#             absorption[wl] = (self.Rp**2 + integral) / (self.Rs**2)


#         count = 0
#         for wl,idx in itertools.product(range(self.nwave),self.iteridx):
#             c,j,k = idx
#             if count < wl:
#                 tau[:] = 0.0 
#                 exptau[:] = 0.0
#                 count += 1    
#             tau[j] += (self.sigma_array[:,wl] * self.X[:,j+k] * self.rho[j+k] * self.dlarray[c])
#             tau[j] += self.Rsig[wl] * rho[j+k] * self.dlarray[c]
#             exptau[j]= exp(-tau[j])
#             integral = 2.0* sum(((self.Rp+self.z)*(1.0-exptau)*self.dz))
#             absorption[wl] = (self.Rp**2 + integral) / (self.Rs**2)
            
#         endtime = time.clock()
#         print 'time ',endtime-starttime