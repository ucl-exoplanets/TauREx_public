'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Transmission class

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

#loading libraries   
import numpy
import copy
import itertools
import time
import ctypes as C
import numpy as np

from library_constants import *
from library_general import *
import logging

import matplotlib.pylab as plt

class transmission():

    def __init__(self, atmosphere, data=None, params=None, ACE=False):

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

        # chemically consistent model
        self.ACE = ACE

        # loading c++ pathintegral library
        self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_transmission.so', mode=C.RTLD_GLOBAL)

        if self.ACE:
            # loading Fortran code for chemically consistent model
            self.ace_lib = C.CDLL('./library/ACE/ACE.so', mode=C.RTLD_GLOBAL)

        # set forward model function
        self.model = self.ctypes_pathintegral

        # set arguments for ctypes libraries
        self.pathintegral_lib.path_integral.argtypes = [C.c_int,
                                                        C.c_int,
                                                        C.c_int,
                                                        C.c_int,
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        C.c_int,
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        C.c_int,
                                                        np.ctypeslib.ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'),
                                                        C.c_int,
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        C.c_int,
                                                        C.c_int,
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                                                        C.c_double,
                                                        C.c_double,
                                                        C.c_void_p]


    def ctypes_pathintegral(self):

        if self.ACE:

            # chemically consistent model
            vector = C.c_double*self.atmosphere.nlayers
            a_apt = vector()
            p_apt = vector()
            t_apt = vector()
            for i in range(nlayers):
               a_apt[i] = self.atmosphere.altitude_profile[i]/1000.
               p_apt[i] = self.atmosphere.pressure_profile[i]/1.e5
               t_apt[i] = self.atmosphere.temperature_profile[i]
            y_out = ((C.c_double * 105) * self.atmosphere.nlayers)()
            self.ace_lib.ACE(C.byref(C.c_int(nlayers)), C.byref(f_x), C.byref(f_x), C.byref(f_x),
                             C.byref(C.c_double(self.atmosphere.He_abund_dex)),
                             C.byref(C.c_double(self.atmosphere.C_abund_dex)),
                             C.byref(C.c_double(self.atmosphere.O_abund_dex)),
                             C.byref(C.c_double(self.atmosphere.N_abund_dex)),
                             C.byref(y_out))
            # y_out has (nlayers, 105). 105 is the total number of molecules computed

            # todo set mixing ratio profiles


        #setting up output array
        absorption = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')


        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_nwngrid,
                                            self.atmosphere.nlayers,
                                            self.atmosphere.nactivegases,
                                            self.atmosphere.ninactivegases,
                                            self.atmosphere.sigma_array_flat,
                                            self.data.sigma_dict['t'],
                                            len(self.data.sigma_dict['t']),
                                            self.atmosphere.sigma_rayleigh_array_flat,
                                            len(self.data.sigma_cia_dict['xsecarr']),
                                            self.atmosphere.cia_idx,
                                            len(self.atmosphere.cia_idx),
                                            self.atmosphere.sigma_cia_array_flat,
                                            self.data.sigma_cia_dict['t'],
                                            len(self.data.sigma_cia_dict['t']),
                                            self.atmosphere.clouds,
                                            self.atmosphere.sigma_clouds_array_flat,
                                            self.atmosphere.clouds_density_profile,
                                            self.atmosphere.density_profile,
                                            self.atmosphere.altitude_profile,
                                            self.atmosphere.active_mixratio_profile.flatten(),
                                            self.atmosphere.inactive_mixratio_profile.flatten(),
                                            self.atmosphere.temperature_profile,
                                            self.atmosphere.planet_radius,
                                            self.params.star_radius,
                                            C.c_void_p(absorption.ctypes.data))

        out = np.zeros((len(absorption)))
        out[:] = absorption

        del(absorption)

        return out





















