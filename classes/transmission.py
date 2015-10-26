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

    ##@profile
    def __init__(self, atmosphere, data=None, params=None, usedatagrid=False):

        logging.info('Initialise object transmission')

        self.__ID__ = 'transmission' #internal class identifier

        if params:
            self.params = params
        else:
            self.params = atmosphere.params #get params object from atmosphere

        if data:
            self.data = data
        else:
            self.data = atmosphere.data    # get data object from atmosphere

        self.atmosphere = atmosphere

        if usedatagrid:
            #use wavelengthgrid of data or internal specgrid defined in data class
            self.lambdagrid = self.data.wavegrid
        else:
            self.lambdagrid = self.data.specgrid
        self.nlambda = len(self.lambdagrid)

        # preload Rayleigh for all gases and the given wavelengths (lambdagrid)
        self.sigma_R_dict = {}
        for gasname in self.params.all_absorbing_gases + self.params.all_inactive_gases:
            if gasname in self.data.sigma_R:
                self.sigma_R_dict[gasname] =  self.data.sigma_R[gasname](self.lambdagrid)

        #calculating collision induced absorption cross sections
        if self.params.in_include_cia:
            self.cia       = self.data.cia
            self.Csig      = self.get_Csig()
        else:
            self.Csig      = zeros((self.nlambda))

        if self.params.in_include_Rayleigh:
            self.Rsig = self.get_Rsig()
        else:
            self.Rsig = zeros((self.nlambda))

        #calculating cloud cross sections
        if self.params.in_include_cld:
            self.Cld_sig = self.get_cloud_sig()
        else:
           self.Cld_sig = zeros((self.nlambda)) # unused but needed to cast to cpp code

        #loading c++ pathintegral library for faster computation
        if self.params.trans_cpp:
            self.cpathlib = C.CDLL('./library/pathintegral.so', mode=C.RTLD_GLOBAL)


        # set forward model function
        self.model = self.cpath_integral
        #self.model = self._path_integral


    #class methods

    def get_dz(self):

        # calculate area of circles of atmosphere
        dz = np.diff(self.atmosphere.z)
        return np.append(dz, dz[-1])

    ##@profile
    def get_Rsig(self):

        # calculating rayleigh scattering cross-sections
        Rsig = zeros((self.nlambda))

        # loop through all gases (absorbing + inactive)
        for gasname in self.params.all_absorbing_gases + self.params.all_inactive_gases:
            if gasname in self.data.sigma_R:
                Rsig += self.atmosphere.get_gas_fraction(gasname) * self.sigma_R_dict[gasname]
        return Rsig

    def get_cloud_sig(self):

        #calculating could cross sections
        a = (128.0* pi**5 * self.atmosphere.clouds_a**6)
        b = (3.0 * self.lambdagrid**4)
        c = ((self.atmosphere.clouds_m**2 -1.0)/(self.atmosphere.clouds_m**2 +2.0))**2

        return a / b * c

    def get_sigma_array(self, temperature):
        # getting sigma array from sigma_dic for given temperature
        return self.data.sigma_dict[find_nearest(self.data.sigma_dict['tempgrid'], temperature)[0]]

    def get_sigma_array_pressure(self):

        # returns  sigma array from sigma_dict for all temperatures, pressures, molecules.

        sigma_dict_arr = np.zeros((len(self.data.sigma_templist),
                                   len(self.data.sigma_preslist),
                                   len(self.params.planet_molec),
                                   self.nlambda))

        for idxtemp, valtemp in enumerate(self.data.sigma_templist):
            for idxpres, valpres in enumerate(self.data.sigma_preslist):
                for idxmol, valmol in enumerate(self.params.planet_molec):
                    sigma_dict_arr[idxtemp, idxpres, idxmol, :] = \
                        self.data.sigma_dict_pres[valtemp][valpres][idxmol]
        return sigma_dict_arr

    def get_sigma_array_c(self):

        #generating 3D sigma_array from sigma_dict for c++ path integral
        tempgrid = self.data.sigma_dict['tempgrid']
        OUT = np.zeros((len(tempgrid), len(self.atmosphere.absorbing_gases),len(self.data.specgrid)), dtype=np.float64)

        c=0
        for t in tempgrid:
            OUT[c,:,:] = self.data.sigma_dict[t]
            c += 1

        return OUT, np.asarray(tempgrid,dtype=np.float64)


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

    # def cget_path_length(self):
    #
    #     nlayers = self.atmosphere.nlayers
    #     czRp = cast2cpp(self.atmosphere.z + self.atmosphere.planet_radius)
    #     dlarray = np.zeros((self.atmosphere.nlayers*self.atmosphere.nlayers), dtype=np.double)
    #     cpath_length = self.cpathlib.cpath_length
    #     cpath_length(C.c_int(nlayers), czRp, C.c_void_p(dlarray.ctypes.data))
    #     return dlarray.reshape(nlayers, nlayers)
    def cget_path_length(self):

        nlayers = self.atmosphere.nlayers
        czRp = cast2cpp(self.atmosphere.z + self.atmosphere.planet_radius)
        dlarray = np.zeros((self.atmosphere.nlayers*self.atmosphere.nlayers), dtype=np.double)
        cpath_length = self.cpathlib.cpath_length
        cpath_length(C.c_int(nlayers), czRp, C.c_void_p(dlarray.ctypes.data))
        return dlarray

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
                cflattened_sigma_arr =  cast2cpp(self.get_sigma_array_pressure().flatten())
                csigma_templist = cast2cpp(self.data.sigma_templist)
                csigma_preslist = cast2cpp(self.data.sigma_preslist)
                cnsigma_templist = C.c_int(len(self.data.sigma_templist))
                cnsigma_preslist = C.c_int(len(self.data.sigma_preslist))
                cpressure_array = cast2cpp(self.atmosphere.P)
                ctemperature = C.c_double(temperature[0])
                
                
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
        clinecount = C.c_int(self.nlambda)
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
        absorption = zeros((self.nlambda),dtype=float64)

        #retrieving function from cpp library
        cpath_int = self.cpathlib.cpath_int

        #running c++ path integral
        cpath_int(csigma_array_2d,csigma_array_3d,cconst_temp,cdlarray,cz,cdz,cRsig,cCsig,cX,crho,cRp,cRs,\
                  clinecount,cnlayers,cn_gas,cn_sig_temp,cInclude_cld,cCld_lowbound,\
                  cCld_upbound,cP_bar,cCld_sig,\
                  cpressure_broadening, cflattened_sigma_arr, csigma_templist, csigma_preslist,\
                  cnsigma_templist, cnsigma_preslist, cpressure_array, ctemperature_array, ctemperature,  \
                  C.c_void_p(absorption.ctypes.data))

        out = zeros((self.nlambda))
        out[:] = absorption

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

    def _path_integral(self, X=None, rho=None, temperature=None, rayleigh=True, absorption=True, clouds=True):

        if X is None:
            X = self.atmosphere.X
        if rho is None:
            rho = self.atmosphere.rho
        if temperature is None:
            temperature = self.params.planet_temp

        nlayers = self.atmosphere.nlayers
        P_bar = self.atmosphere.P_bar
        dlarray, iteridx = self._get_path_length()
        dz = self.get_dz()

        # calculating rayleigh scattering cross sections
        if self.params.in_include_Rayleigh:
            self.Rsig = self.get_Rsig()
        else:
            self.Rsig = zeros((self.nlambda))

        #selecting correct sigma_array for temperature
        sigma_array = self.get_sigma_array(temperature)

        #setting up arrays
        absorption = zeros((self.nlambda))
        tau = zeros((nlayers))
        Rtau = zeros((nlayers))
        Ctau = zeros((nlayers))
        cld_tau = zeros((nlayers))
        exptau = zeros((nlayers))

        molnum = len(X[:,0])

        for wl in range(self.nlambda): # loop through wavelenghts

            tau[:] = 0.0
            exptau[:] = 0.0
            Rtau[:] = 0.0
            cld_tau[:] = 0.0

            for c, j, k in iteridx:

                # optical depth due to gasses
                for i in range(molnum):
                    tau[j] += (sigma_array[i,wl] * X[i,j+k] * rho[j+k] * dlarray[c])

                #optical depth due to Rayleigh scattering
                Rtau[j] = self.Rsig[wl] * rho[j+k] * dlarray[c]

                # calculating CIA optical depth
                Ctau[j] = self.Csig[wl] * (rho[j+k] **2) * dlarray[c]

                #Calculating cloud opacities, single cloud layer implementation only. can be upgraded or ... not
                if self.params.in_include_cld:
                    if P_bar[k+j] < self.atmosphere.clouds_upper_P and P_bar[k+j] > self.atmosphere.clouds_lower_P:
                        # = log(cloud density), assuming linear decrease with decreasing log pressure
                        # following Ackerman & Marley (2001), Fig. 6
                        cld_log_rho = interp_value(np.log(P_bar[k+j]),np.log(self.atmosphere.clouds_lower_P),np.log(self.atmosphere.clouds_upper_P),-6.0,-1.0)
                        cld_tau[j] += ( self.Cld_sig[wl] * (dlarray[c]*1.0e2) * (exp(cld_log_rho)*1.0e-6) )    # convert path lenth from m to cm, and density from g m^-3 to g cm^-3

                #adding all taus together
                tau[j] += Rtau[j]
                tau[j] += Ctau[j]
                tau[j] += cld_tau[j]

                exptau[j]= exp(-tau[j])

            integral = 2.0* sum(((self.atmosphere.planet_radius+self.atmosphere.z)*(1.0-exptau)*dz))
            absorption[wl] = (self.atmosphere.planet_radius**2 + integral) / (self.params.star_radius**2)

        return absorption

