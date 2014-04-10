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
import numpy
from numpy import *
from pylab import *
import itertools
import time
import ctypes as C
from library.library_general import *

class transmission(object):

#initialisation
    def __init__(self,params,data):
        #loading data
        self.params        = params
        self.z             = data['pta'][:,2]
        self.Rp            = params['planet_radius']
        self.Rs            = params['star_radius']
        self.nlayers       = data['nlayers']
        self.rho           = data['rho']
        self.n_gas         = data['ngas']
        self.sigma_array   = data['sigma_array']
        self.wavegrid      = data['wavegrid']
        self.nwave         = data['nwave']
        self.nspecgrid     = data['nspecgrid']
        self.X             = data['X']
        self.atmosphere    = data['atmosphere']

        #calculating optical path lengths
        self.dlarray,self.iteridx = self.get_path_length()
        self.dz            = self.get_dz()
 
        #calculating rayleigh scattering cross sections
        self.mol_list      = self.atmosphere['mol'].keys()
        self.ray_thres     = 50.0 * (self.atmosphere['mol'][self.mol_list[0]]['radius'] * 1.0e6);
        self.Rsig          = self.get_Rsig()

        #calculating collision induced absorption cross sections
        if params.in_include_cia == True:
            self.cia       = data['cia']
            self.Csig      = self.get_Csig()
        else:
            self.Csig      = zeros((self.nwave))
        
        #loading c++ pathintegral library for faster computation
        if params.trans_cpp == True:
            self.cpathlib = C.cdll.LoadLibrary('./library/pathintegral.so')

     
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


#class methods

    def get_path_length(self):
        #calculates the layerlength 
        dlarray = []
        jl = []
        kl = []
        cl = []
        count = 0
        for j in range(self.nlayers-1):    
            for k in range(1,self.nlayers-j):
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
        dz[-1] = dz[-2]
        return dz
    
    def get_Rsig(self):
    #calculating rayleigh scattering cross-sections
        Rsig = zeros((self.nwave))
        count = 0
        for wl in self.wavegrid:
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



    def path_integral(self, X = None, rho = None):

        if X == None:
            X = self.X
        if rho == None:
            rho = self.rho
        
        #setting up arrays
        absorption = zeros((self.nspecgrid))
        tau = zeros((self.nlayers))  
        Rtau = zeros((self.nlayers))
        Ctau = zeros((self.nlayers))
        exptau = zeros((self.nlayers))

        molnum = len(X[:,0])
        #running loop over wavelengths     
        for wl in range(self.nspecgrid):
            tau[:] = 0.0 
            exptau[:] = 0.0
            Rtau[:] = 0.0
         
            for c,j,k in self.iteridx:
                #optical depth due to gasses
                for i in range(molnum):
                    tau[j] += (self.sigma_array[i,wl] * X[i,j+k] * rho[j+k] * self.dlarray[c])
                
                Rtau[j] = self.Rsig[wl] * rho[j+k] * self.dlarray[c] #optical depth due to Rayleigh scattering
                Ctau[j] = self.Csig[wl] * (rho[j+k] **2) * self.dlarray[c]  # calculating CIA optical depth
                #adding all taus together
                tau[j] += Rtau[j]
                tau[j] += Ctau[j]
 
                exptau[j]= exp(-tau[j])
  
            integral = 2.0* sum(((self.Rp+self.z)*(1.0-exptau)*self.dz))
#             print self.z, self.dz, exptau
            
            absorption[wl] = (self.Rp**2 + integral) / (self.Rs**2)

        return absorption
        
        
    def cpath_integral(self, X = None, rho = None):
        if X == None:
            X = self.X
        if rho == None:
            rho = self.rho        

        #casting changing arrays to c++ pointers
        Xs1,Xs2 = shape(X)
        Xnew = zeros((Xs1+1,Xs2))
        Xnew[:-1,:] = X
        cX, cs1,cs2 = cast2cpp(Xnew)
        crho, cs1 = cast2cpp(rho)
        
        #casting fixed arrays and variables to c++ pointers
        csigma_array, cs1,cs2 = cast2cpp(self.sigma_array)
        cdlarray, cs1 = cast2cpp(self.dlarray)
        znew = zeros((len(self.z)))
        znew[:] = self.z
        cz, cs1   = cast2cpp(znew)
        cdz, cs1  = cast2cpp(self.dz) 
        cRsig, cs1 = cast2cpp(self.Rsig)
        cCsig, cs1 = cast2cpp(self.Csig)
        cRp = C.c_double(self.Rp)
        cRs = C.c_double(self.Rs)
        clinecount = C.c_int(self.nspecgrid)
        cnlayers = C.c_int(self.nlayers)
        cn_gas = C.c_int(len(X[:,0]))
        
        #setting up output array
        absorption = zeros((self.nspecgrid),dtype=float64)
        
        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int
        
        #running c++ path integral
        cpath_int(csigma_array,cdlarray,cz,cdz,cRsig,cCsig,cX,crho,cRp,cRs,\
                  clinecount,cnlayers,cn_gas,C.c_void_p(absorption.ctypes.data))

        return absorption   
        
        
    
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
        wl=wl *1.0e-6               #convert wavelengths to m
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
# #                     for l in range(self.n_gas):
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