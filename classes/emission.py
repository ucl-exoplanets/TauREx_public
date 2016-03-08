'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Emission class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries   
import copy
import itertools
import time
import ctypes as C
import numpy as np

from library_constants import *
from library_general import *
import logging

import matplotlib.pylab as plt

class emission(object):

    def __init__(self, atmosphere, data=None, params=None):

        logging.info('Initialise object emission')

        if params:
            self.params = params
        else:
            self.params = atmosphere.params #get params object from atmosphere
        if data:
            self.data = data
        else:
            self.data = atmosphere.data    # get data object from atmosphere

        self.atmosphere = atmosphere

        #loading c++ pathintegral library for faster computation
        self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_emission.so', mode=C.RTLD_GLOBAL)

        # set forward model function
        self.model = self.ctypes_pathintegral

        #retrieving function from cpp library
        self.pathintegral_lib.path_integral.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         C.c_int,
                                                         C.c_int,
                                                         C.c_int,
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         C.c_int,
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         C.c_double,
                                                         C.c_double,
                                                         np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                         C.c_void_p]


    def ctypes_pathintegral(self, return_tau=False):

        #setting up output array
        FpFs = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_wngrid,
                                             self.atmosphere.int_nwngrid,
                                             self.atmosphere.nlayers,
                                             self.atmosphere.nactivegases,
                                             self.atmosphere.sigma_array_flat,
                                             self.data.sigma_dict['t'],
                                             len(self.data.sigma_dict['t']),
                                             self.atmosphere.density_profile,
                                             self.atmosphere.altitude_profile,
                                             self.atmosphere.active_mixratio_profile.flatten(),
                                             self.atmosphere.temperature_profile,
                                             self.atmosphere.planet_radius,
                                             self.params.star_radius,
                                             self.atmosphere.data.star_sed,
                                             C.c_void_p(FpFs.ctypes.data),
                                             C.c_void_p(tau.ctypes.data))



        out = np.zeros((len(FpFs)))
        out[:] = FpFs

        if return_tau:

            tauout = np.zeros((self.atmosphere.int_nwngrid, self.atmosphere.nlayers))
            count = 0
            for i in range(self.atmosphere.int_nwngrid):
                for j in range(1, self.atmosphere.nlayers):
                    tauout[i,j] = tau[count]
                    count += 1
            tauout = np.fliplr(np.rot90(tauout))

            del(FpFs)
            del(tau)

            return tauout

        else:

            out = np.zeros((len(FpFs)))
            out[:] = FpFs

            del(FpFs)
            del(tau)

            return out