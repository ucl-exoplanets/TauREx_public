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

    def __init__(self, atmosphere):

        logging.info('Initialise object transmission')


        self.params = atmosphere.params #get params object from atmosphere
        self.data = atmosphere.data    # get data object from atmosphere
        self.atmosphere = atmosphere

        if self.params.gen_ace:
            # loading Fortran code for chemically consistent model
            self.ace_lib = C.CDLL('./library/ACE/ACE.so', mode=C.RTLD_GLOBAL)

        if self.params.in_opacity_method in ['xsec_lowres', 'xsec_highres', 'xsec_sampled', 'xsec']: # using cross sections

            # loading c++ pathintegral library
            if self.atmosphere.nthreads > 1:
                # load openmp version. Remember to set OMP_NUM_THREADS to number of threads
                logging.info('Load openmp version of transmission model')
                self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_transmission_parallel_xsec.so', mode=C.RTLD_GLOBAL)
            else:
                logging.info('Load single-core version of transmission model')
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
                C.c_int, # params.atm_mie
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
                np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                C.c_double,
                C.c_double,
                C.c_void_p,
                C.c_void_p]


    def ctypes_pathintegral_xsec(self, return_tau=False, mixratio_mask=False):

        if self.params.gen_ace:
            self.atmosphere.set_ACE(mixratio_mask)

        #setting up output array
        absorption = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_nwngrid,
                                            self.atmosphere.nlayers,
                                            self.atmosphere.nactivegases,
                                            self.atmosphere.ninactivegases,
                                            self.params.atm_rayleigh,
                                            self.params.atm_mie,
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
                                            self.atmosphere.clouds_pressure,
                                            self.atmosphere.sigma_mie_array,
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


    def ctypes_pathintegral_ktab(self, return_tau=False, mixratio_mask=False):

        if self.params.gen_ace:
            self.atmosphere.set_ACE(mixratio_mask)


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
                                            self.atmosphere.clouds_pressure,
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