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
import copy
import itertools
import time
import ctypes as C
import numpy as np

from library_constants import *
from library_general import *
from library_transmission import *
import logging

import matplotlib.pylab as plt

class transmission():

    def __init__(self, atmosphere, data=None, params=None):

        logging.info('Initialise object transmission')

        if params:
            self.params = params
        else:
            self.params = atmosphere.params #get params object from atmosphere
        if data:
            self.data = data
        else:
            self.data = atmosphere.data    # get data object from atmosphere

        self.atmosphere = atmosphere

        # todo IMPROVE RAYLEIGH SCATTERING
        # preload Rayleigh for all gases and the given wavelengths (lambdagrid)
        self.sigma_R_dict = {}
        for gasname in self.params.all_absorbing_gases + self.params.all_inactive_gases:
            if gasname in self.data.sigma_R:
                self.sigma_R_dict[gasname] =  self.data.sigma_R[gasname](self.data.int_wngrid)
        if self.params.in_include_Rayleigh:
            self.Rsig = self.get_Rsig()
        else:
            self.Rsig = zeros((self.data.int_nwngrid))

        #calculating collision induced absorption cross sections
        if self.params.in_include_cia:
            self.cia = self.data.cia
            self.Csig = self.get_Csig()
        else:
            self.Csig = zeros((self.data.int_nwngrid))

        #calculating cloud cross sections
        if self.params.in_include_cld:
            self.Cld_sig = self.get_cloud_sig()
        else:
           self.Cld_sig = zeros((self.data.int_nwngrid)) # unused but needed to cast to cpp code

        #loading c++ pathintegral library for faster computation
        if self.params.trans_cpp:
            self.cpathlib = C.CDLL('./library/cpath_pathintegral.so', mode=C.RTLD_GLOBAL)

        # set forward model function
        self.cpath_load_vars() # preload variables in memory for cpath
        self.model = self.cpath_pathintegral

    #class methods

    def get_dz(self):
        # calculate area of circles of atmosphere
        dz = np.diff(self.atmosphere.z)
        return np.append(dz, dz[-1])

    def get_Rsig(self):

        # calculating rayleigh scattering cross-sections
        Rsig = zeros((self.data.int_nwngrid))
        # loop through all gases (absorbing + inactive)
        for gasname in self.params.all_absorbing_gases + self.params.all_inactive_gases:
            if gasname in self.data.sigma_R:
                Rsig += self.atmosphere.get_gas_fraction(gasname) * self.sigma_R_dict[gasname]
        return Rsig

    def get_cloud_sig(self):
        #calculating could cross sections
        a = (128.0* pi**5 * self.atmosphere.clouds_a**6)
        b = (3.0 * (10000./self.data.obs_wngrid)**4)
        c = ((self.atmosphere.clouds_m**2 -1.0)/(self.atmosphere.clouds_m**2 + 2.0))**2
        return a / b * c

    def get_Csig(self):
        '''

        calculating collisional induced absorption

        Optical depth given by tau = alpha * L * c_1 * c_2 ; alpha = abs coeff (cm^5 mol^-2), L = path length (cm), c_i = concentration of collider i (mol cm^-3)
        IN: CIA coeffs in (cm^-1 amagat^-2), grid wavelength (in um), path length (in m) and total number density dz (in m^-3)
        OUT: H2-H2 collision-induced absorption coefficient (in m^5 mol^-2)
        '''
    #     AMA=2.68676e+25 #Amagat (molecules m^-3)
    #     conv_factor = 1.0 / pow((AMA*1.0e-6),2) #converting from (cm^-1 amagat^-2) to (cm^5 mol^-2)...
    #     alpha = coeff * conv_factor
    #     alpha *= (amount**2) * 1.0e-10 #e.g. composition 85% H2, and convert from cm^5 to m^5

        sigma = self.cia[:,1] * (self.atmosphere.get_gas_fraction('H2')**2) * 1.385294e-49
        return sigma


    def cpath_path_length(self):

        zRp = self.atmosphere.z + self.atmosphere.planet_radius
        czRp = cast2cpp(zRp)
        dlarray = np.zeros((self.atmosphere.nlayers*self.atmosphere.nlayers), dtype=np.double)
        cpath_length = self.cpathlib.path_length
        cpath_length(C.c_int(self.atmosphere.nlayers), czRp, C.c_void_p(dlarray.ctypes.data))

        return dlarray

    def cpath_load_vars(self):



        self.cpath_nwngrid = C.c_int(self.data.int_nwngrid)
        self.cpath_nlayers = C.c_int(self.atmosphere.nlayers)
        self.cpath_nmol = C.c_int(self.data.sigma_nmol)
        self.cpath_sigma_array = cast2cpp(self.data.sigma_array)
        self.cpath_sigma_np = cast2cpp(self.data.sigma_np.flatten())
        self.cpath_sigma_nt = cast2cpp(self.data.sigma_nt.flatten())
        self.cpath_sigma_p = cast2cpp(self.data.sigma_p.flatten())
        self.cpath_sigma_t = cast2cpp(self.data.sigma_t.flatten())
        self.cpath_dlarray = cast2cpp(self.cpath_path_length())
        self.cpath_z = cast2cpp(self.get_dz())
        self.cpath_dz = cast2cpp(self.atmosphere.z)
        #self.cpath_rho = cast2cpp(self.atmosphere.rho)
        #self.cpath_X = cast2cpp(self.atmosphere.X)
        #self.cpath_T = cast2cpp(self.atmosphere.T)
        #self.cpath_P = cast2cpp(self.atmosphere.P)
        #self.cpath_Rp = cast2cpp(self.atmosphere.planet_radius)
        self.cpath_Rs = cast2cpp(self.params.star_radius)
    
    def cpath_update_vars(self):
        # todo should update depending on fitted params
        self.cpath_X = cast2cpp(self.atmosphere.X)
        self.cpath_T = cast2cpp(self.atmosphere.T)
        self.cpath_P = cast2cpp(self.atmosphere.P)
        self.cpath_Rp = cast2cpp(self.atmosphere.planet_radius)
        self.cpath_rho = cast2cpp(self.atmosphere.rho)
    
    def cpath_pathintegral(self):
        
        self.cpath_update_vars()
    
        #setting up output array
        absorption = zeros((self.data.int_nwngrid),dtype=float64)
        
        #retrieving function from cpp library
        cpath_integral = self.cpathlib.path_integral
        
        #running c++ path integral
        cpath_integral(self.cpath_nwngrid,
                       self.cpath_nlayers,
                       self.cpath_nmol,
                       self.cpath_sigma_array,
                       self.cpath_sigma_np,
                       self.cpath_sigma_nt,
                       self.cpath_sigma_p,
                       self.cpath_sigma_t,
                       self.cpath_dlarray,
                       self.cpath_z,
                       self.cpath_dz,
                       self.cpath_rho,
                       self.cpath_X,
                       self.cpath_T,
                       self.cpath_P,
                       self.cpath_Rp,
                       self.cpath_Rs,
                       C.c_void_p(absorption.ctypes.data))

        print " there are ", self.atmosphere.nlayers

        print self.atmosphere.T

        print 'end'
        exit()

    #@profile
    def cpath_integral(self, X=None, rho=None, temperature=None):

        if X is None:
            X = self.atmosphere.X
        if rho is None:
            rho = self.atmosphere.rho
        if temperature is None:
            temperature = self.atmosphere.T

        dz = self.get_dz()
        dlarray = self.cget_path_length()

        #selecting correct sigma_array for temperature array
        cpressure_broadening = C.c_int(0)

        cflattened_sigma_arr =  cast2cpp(np.zeros(0))
        csigma_templist = cast2cpp(np.zeros(0))
        csigma_preslist = cast2cpp(np.zeros(0))
        cnsigma_templist = C.c_int(0)
        cnsigma_preslist = C.c_int(0)
        cpressure_array = cast2cpp(np.zeros(0))

        ctemperature = C.c_double(temperature[0])
        ctemperature_array = cast2cpp(np.zeros(0))

        sigma_array_3d = np.zeros(0)
        sigma_array_2d = np.zeros(0)
        cn_sig_temp= C.c_int()
        cconst_temp = C.c_int(0)

        if len(np.unique(temperature)) == 1:

            # constant T with altitude
            cconst_temp = C.c_int(1)
            # cast first value of T profile (it's isothermal..)
            ctemperature = C.c_double(temperature[0])

            if self.params.in_use_P_broadening:
                
                cpressure_broadening = C.c_int(1)
                flattened_sigma_arr = self.get_sigma_array_pressure().flatten()
                cflattened_sigma_arr =  (C.c_double * len(flattened_sigma_arr))(*flattened_sigma_arr)



                csigma_templist = cast2cpp(self.data.sigma_templist)
                csigma_preslist = cast2cpp(self.data.sigma_preslist)
                cnsigma_templist = C.c_int(len(self.data.sigma_templist))
                cnsigma_preslist = C.c_int(len(self.data.sigma_preslist))
                cpressure_array = cast2cpp(self.atmosphere.P)
                ctemperature = C.c_double(temperature[0])
                cn_sig_temp = cnsigma_templist # again?

            else:
                sigma_array_2d = self.get_sigma_array(temperature[0])
        else:

            logging.debug('Variable T with altitude. Get 3d sigma array and set other variables')

            # variable T with altitude. Get 3d sigma array and set other variables
            # note this does not include pressure broadening
            sigma_array_3d, tempgrid = self.get_sigma_array_c()

            sigma_array_3d = sigma_array_3d.flatten()


            ctemperature_array = (C.c_double * len(self.atmosphere.T.tolist()))(*self.atmosphere.T.tolist())
            #ctemperature_array = cast2cpp(self.atmosphere.T) # crazy stuff, here cast2cpp doesn't work !!!
            csigma_templist = cast2cpp(tempgrid) # temperature grid of cross sections


            cnsigma_templist= C.c_int(len(tempgrid))
            cn_sig_temp= C.c_int(len(tempgrid)) # legacy???

        #casting to c++ pointers
        csigma_array_2d = cast2cpp(sigma_array_2d)
        csigma_array_3d = cast2cpp(sigma_array_3d)

        Xs1, Xs2 = shape(X)
        Xnew = zeros((Xs1+1, Xs2))
        Xnew[:-1,:] = X
        cX = cast2cpp(Xnew)
        crho = cast2cpp(rho)

        cdlarray = cast2cpp(dlarray)
        znew = zeros((len(self.atmosphere.z)))
        znew[:] = self.atmosphere.z
        cz = cast2cpp(znew)
        cdz = cast2cpp(dz)
        cRsig = cast2cpp(self.Rsig)

        cCsig = cast2cpp(self.Csig)

        cRp = C.c_double(self.atmosphere.planet_radius)
        cRs = C.c_double(self.params.star_radius)
        clinecount = C.c_int(self.data.int_nwngrid)
        cnlayers = C.c_int(self.atmosphere.nlayers)

        cn_gas = C.c_int(len(X[:,0]))

        #setting and casting cloud paramters
        cP_bar = cast2cpp(self.atmosphere.P_bar)
        cCld_sig = cast2cpp(self.Cld_sig)
        if self.params.in_include_cld:
            cInclude_cld = C.c_int(1)
            cCld_lowbound = C.c_double(self.atmosphere.clouds_lower_P)
            cCld_upbound = C.c_double(self.atmosphere.clouds_upper_P)

        else:
            cInclude_cld = C.c_int(0)
            cCld_lowbound = C.c_double(0)
            cCld_upbound = C.c_double(0)

        #setting up output array
        absorption = zeros((self.data.int_nwngrid),dtype=float64)

        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int

        #running c++ path integral
        cpath_int(cflattened_sigma_arr, csigma_array_2d,csigma_array_3d,cconst_temp,cdlarray,cz,cdz,cRsig,cCsig,cX,crho,cRp,cRs,\
                  clinecount,cnlayers,cn_gas,cn_sig_temp,cInclude_cld,cCld_lowbound,\
                  cCld_upbound,cP_bar,cCld_sig,\
                  cpressure_broadening, csigma_templist, csigma_preslist,\
                  cnsigma_templist, cnsigma_preslist, cpressure_array, ctemperature_array, ctemperature,  \
                  C.c_void_p(absorption.ctypes.data))

        out = zeros((self.data.int_nwngrid))
        out[:] = absorption

        del(dlarray)
        del(absorption)


        return out

    # deprecated functions

    def _get_path_length(self):

        nlayers = self.atmosphere.nlayers
        z = self.atmosphere.z
        Rp = self.atmosphere.planet_radius

        dlarray = []
        jl = []
        kl = []
        cl = []
        count = 0
        for j in range(nlayers): # loop through atmosphere layers
            for k in range(1,nlayers-j): # loop through each layer to sum up path length
                dl = 2.0 * (sqrt(pow((Rp + z[k+j]),2) - pow((Rp + z[j]),2)) -
                            sqrt(pow((Rp + z[k-1+j]),2) - pow((Rp + z[j]),2)))
                dlarray.append(dl)
                jl.append(j)
                kl.append(k)
                cl.append(count)
                count += 1

        iteridx = zip(cl,jl,kl)
        dlarray = asarray(dlarray, dtype=np.float)
        return dlarray, iteridx


