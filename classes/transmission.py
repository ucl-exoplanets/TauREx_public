'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Transmission class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries   
import numpy
import ctypes as C
import numpy as np

from library_constants import *
from library_general import *
import logging

try:
    import cPickle as pickle
except:
    import pickle

import matplotlib.pylab as plt

class transmission(object):

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

        if self.params.in_opacity_method in ['xsec_highres', 'xsec_lowres']: # using cross sections

            # loading c++ pathintegral library
            if self.atmosphere.nthreads > 1:
                # load openmp version. Remember to set OMP_NUM_THREADS to number of threads
                self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_transmission_parallel_xsec.so', mode=C.RTLD_GLOBAL)
            else:
                self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_transmission_xsec.so', mode=C.RTLD_GLOBAL)

            self.model = self.ctypes_pathintegral_xsec
            # set arguments for ctypes libraries
            self.pathintegral_lib.path_integral.argtypes = [
                C.c_int,
                C.c_int,
                C.c_int,
                C.c_int,
                C.c_int,
                C.c_int,
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                C.c_double,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_double,
                C.c_double,
                C.c_void_p,
                C.c_void_p]


        elif self.params.in_opacity_method in ['ktab', 'ktable', 'ktables']: # using k tables
            # loading c++ pathintegral library
            self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_transmission_ktab.so', mode=C.RTLD_GLOBAL)
            self.model = self.ctypes_pathintegral_ktab
            # set arguments for ctypes libraries
            self.pathintegral_lib.path_integral.argtypes = [
                C.c_int, # atmosphere.int_nwngrid
                C.c_int, # atmosphere.nlayers
                C.c_int, # atmosphere.nactivegases
                C.c_int, # atmosphere.ninactivegases
                C.c_int, # params.atm_rayleigh
                C.c_int, # params.atm_cia
                C.c_int, # params.atm_clouds
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # atmosphere.ktables_array_flat
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # data.ktable_dict['t']
                C.c_int, # len(data.ktable_dict['t'])
                C.c_int, # data.ktable_dict['ngauss']
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # data.ktable_dict['weights']
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_int,
                C.c_double,
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_double,
                C.c_double,
                C.c_void_p,
                C.c_void_p]




        # set forward model function
        #self.model = self.ctypes_pathintegral

        #self.model = self.linebyline_pathintegral


    def ktables_pathintegral(self):

        z = self.atmosphere.altitude_profile
        dz = np.diff(z)
        dz = np.append(dz, dz[-1])
        nlayers = self.atmosphere.nlayers

        rp = self.atmosphere.planet_radius
        density = self.atmosphere.density_profile
        mixratio = self.atmosphere.active_mixratio_profile
        dlarray = np.zeros((nlayers, nlayers))

        for j in range(nlayers):
            for k in range(1, nlayers - j):
                p = np.power(z[j] + rp, 2)
                dlarray[j, k] = 2.*np.sqrt(np.power(z[k+j]+rp, 2) - p) - np.sqrt(np.power(z[k-1+j]+rp, 2) - p)

        # load k table
        #ktable = []
        #ktable.append(pickle.load(open('/Users/marco/Dropbox/repos/TauREx/Output/ktables/ktable_h2o_0.001bar_1500K.kcoeff.pickle')))
        #ktable.append(pickle.load(open('/Users/marco/Dropbox/repos/TauREx/Output/ktables/ktable_co_1bar_1200K.kcoeff.pickle')))

        mixratio[1, :] = 1e-3

        bins = ktable[0][:,0]

        ngauss = 20
        gauss = np.polynomial.legendre.leggauss(ngauss)
        weights = gauss[1]

        absorption = np.zeros(len(bins))

        # for bin_idx, bin_val in enumerate(bins):
        #     integral = 0
        #     for j in range(nlayers):
        #         transfull = 1.
        #         for gasid in range(len(ktable)):
        #             kcoeff = ktable[gasid][:,1:]
        #             transtmp = 0
        #             for kc_idx in range(ngauss):
        #                 tautmp = 0
        #                 for k in range(1, nlayers-j):
        #                     tautmp += kcoeff[bin_idx, kc_idx] * mixratio[gasid, j+k] * density[j+k] * dlarray[j, k]
        #                 transtmp += np.exp(- tautmp ) * weights[kc_idx]/2.
        #             transfull *= transtmp
        #         integral += (rp + z[j])*(1-transfull)*dz[j]
        #     integral *= 2.
        #     absorption[bin_idx] = (rp**2 + integral)/(self.params.star_radius**2)

        for bin_idx, bin_val in enumerate(bins):
            integral = 0
            for j in range(nlayers):
                transfull = 1.
                for gasid in range(len(ktable)):
                    kcoeff = ktable[gasid][:,1:]
                    transtmp = 0
                    for kc_idx in range(ngauss):
                        tautmp = 0
                        for k in range(1, nlayers-j):
                            tautmp += kcoeff[bin_idx, kc_idx] * mixratio[gasid, j+k] * density[j+k] * dlarray[j, k]
                        transtmp += np.exp(- tautmp ) * weights[kc_idx]/2.


                    transfull *= transtmp
                integral += (rp + z[j])*(1-transfull)*dz[j]
            integral *= 2.
            absorption[bin_idx] = (rp**2 + integral)/(self.params.star_radius**2)


        out = np.zeros((len(absorption), 2))
        out[:,1] = absorption
        out[:,0] = 10000./bins

        np.savetxt('Output/ktables/spectrum_kcoeff_0.001bar.dat', out)
        exit()

    def linebyline_pathintegral(self):

        z = self.atmosphere.altitude_profile
        dz = np.diff(z)
        dz = np.append(dz, dz[-1])
        nlayers = self.atmosphere.nlayers

        rp = self.atmosphere.planet_radius
        density = self.atmosphere.density_profile
        mixratio = self.atmosphere.active_mixratio_profile
        dlarray = np.zeros((nlayers, nlayers))

        # mixratio[1, :] = 1e-3

        for j in range(nlayers):
            for k in range(1, nlayers - j):
                p = np.power(z[j] + rp, 2)
                dlarray[j, k] = 2.*np.sqrt(np.power(z[k+j]+rp, 2) - p) - np.sqrt(np.power(z[k-1+j]+rp, 2) - p)


        #load xsec
        xsec = []
        #xsec.append(pickle.load(open('/Users/marco/Dropbox/repos/TauREx/Output/ktables/1H2-16O__0-30000_1500_0.001_0.01.sig.pickle')))
        #xsec.append(pickle.load(open('/Users/marco/Dropbox/repos/TauREx/Output/ktables/12C-16O__0-30000_1200_1.0_0.01.sig.pickle')))

        absorption = np.zeros(len(xsec[0][:,0]))

        for bin_idx, bin_val in enumerate(xsec[0][:,0]):
            integral = 0
            for j in range(nlayers):
                tautmp = 0
                for k in range(1, nlayers-j):
                    for gasid in range(len(xsec)):
                        tautmp += xsec[gasid][bin_idx, 1] * mixratio[gasid, j+k] * density[j+k] * dlarray[j, k]

                integral += (rp + z[j])*(1-np.exp(- tautmp))*dz[j]

            integral *= 2.
            absorption[bin_idx] = (rp**2 + integral)/(self.params.star_radius**2)


        out = np.zeros((len(absorption), 2))
        out[:,1] = absorption
        out[:,0] = 10000./xsec[0][:,0]

        np.savetxt('Output/ktables/spectrum_linebyline_0.001bar.dat', out)
        exit()



    def ctypes_pathintegral_xsec(self, return_tau=False, mixratio_mask=False):

        if self.params.gen_ace:
            self.set_ACE(mixratio_mask)

        #setting up output array
        absorption = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_nwngrid,
                                            self.atmosphere.nlayers,
                                            self.atmosphere.nactivegases,
                                            self.atmosphere.ninactivegases,
                                            self.params.atm_rayleigh,
                                            self.params.atm_cia,
                                            self.params.atm_clouds,
                                            self.atmosphere.sigma_array_flat,
                                            self.data.sigma_dict['t'],
                                            len(self.data.sigma_dict['t']),
                                            self.atmosphere.sigma_rayleigh_array_flat,
                                            len(self.data.sigma_cia_dict['xsecarr']),
                                            np.asarray(self.atmosphere.cia_idx, dtype=np.float),
                                            len(self.atmosphere.cia_idx),
                                            self.atmosphere.sigma_cia_array_flat,
                                            self.data.sigma_cia_dict['t'],
                                            len(self.data.sigma_cia_dict['t']),
                                            self.atmosphere.clouds_topP,
                                            self.atmosphere.pressure_profile,
                                            self.atmosphere.density_profile,
                                            self.atmosphere.altitude_profile,
                                            self.atmosphere.active_mixratio_profile.flatten(),
                                            self.atmosphere.inactive_mixratio_profile.flatten(),
                                            self.atmosphere.temperature_profile,
                                            self.atmosphere.planet_radius,
                                            self.params.star_radius,
                                            C.c_void_p(absorption.ctypes.data),
                                            C.c_void_p(tau.ctypes.data))

        out = np.zeros((len(absorption)))
        out[:] = absorption

        if return_tau:

            tauout = np.zeros((self.atmosphere.int_nwngrid, self.atmosphere.nlayers))
            count = 0
            for i in range(self.atmosphere.int_nwngrid):
                for j in range(self.atmosphere.nlayers):
                    tauout[i,j] = tau[count]
                    count += 1
            tauout = np.fliplr(np.rot90(tauout))

            del(absorption)
            del(tau)

            return tauout

        else:

            del(absorption)
            del(tau)

            return out


    def ctypes_pathintegral_ktab(self, return_tau=False, mixratio_mask=False):

        if self.params.gen_ace:
            self.set_ACE(mixratio_mask)

        #setting up output array
        absorption = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_nwngrid,
                                            self.atmosphere.nlayers,
                                            self.atmosphere.nactivegases,
                                            self.atmosphere.ninactivegases,
                                            self.params.atm_rayleigh,
                                            self.params.atm_cia,
                                            self.params.atm_clouds,
                                            self.atmosphere.ktables_array_flat,
                                            self.data.ktable_dict['t'],
                                            len(self.data.ktable_dict['t']),
                                            self.data.ktable_dict['ngauss'],
                                            self.data.ktable_dict['weights'],
                                            self.atmosphere.sigma_rayleigh_array_flat,
                                            len(self.data.sigma_cia_dict['xsecarr']),
                                            np.asarray(self.atmosphere.cia_idx, dtype=np.float),
                                            len(self.atmosphere.cia_idx),
                                            self.atmosphere.sigma_cia_array_flat,
                                            self.data.sigma_cia_dict['t'],
                                            len(self.data.sigma_cia_dict['t']),
                                            self.atmosphere.clouds_topP,
                                            self.atmosphere.pressure_profile,
                                            self.atmosphere.density_profile,
                                            self.atmosphere.altitude_profile,
                                            self.atmosphere.active_mixratio_profile.flatten(),
                                            self.atmosphere.inactive_mixratio_profile.flatten(),
                                            self.atmosphere.temperature_profile,
                                            self.atmosphere.planet_radius,
                                            self.params.star_radius,
                                            C.c_void_p(absorption.ctypes.data),
                                            C.c_void_p(tau.ctypes.data))

        out = np.zeros((len(absorption)))
        out[:] = absorption

        if return_tau:

            tauout = np.zeros((self.atmosphere.int_nwngrid, self.atmosphere.nlayers))
            for i in range(self.atmosphere.int_nwngrid):
                for j in range(self.atmosphere.nlayers):
                    tauout[i,j] = tau[i + j*self.atmosphere.int_nwngrid]

            tauout = np.fliplr(np.rot90(tauout))

            del(absorption)
            del(tau)

            return tauout

        else:

            del(absorption)
            del(tau)

            return out



