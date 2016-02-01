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

        #loading c++ pathintegral library for faster computation
        self.cpathlib = C.CDLL('./library/ctypes_pathintegral.so', mode=C.RTLD_GLOBAL)
        C._reset_cache()

        # preload variables in memory for cpath
        self.ctypes_load_vars()

        # set forward model function
        self.model = self.ctypes_pathintegral

    # class methods

    def ctypes_load_vars(self):

        # load variables that won't change during fitting
        self.ctypes_nwngrid = C.c_int(self.atmosphere.int_nwngrid)
        self.ctypes_nlayers = C.c_int(self.atmosphere.nlayers)
        self.ctypes_nactive = C.c_int(self.atmosphere.nactivegases)
        self.ctypes_ninactive = C.c_int(self.atmosphere.ninactivegases)
        self.ctypes_sigma_array = cast2cpp(self.atmosphere.sigma_array_flat)
        self.ctypes_sigma_temp = cast2cpp(self.data.sigma_dict['t'])
        self.ctypes_sigma_ntemp = C.c_int(len(self.data.sigma_dict['t']))
        self.ctypes_sigma_cia_array = cast2cpp(self.atmosphere.sigma_cia_array_flat)
        self.ctypes_cia_npairs = C.c_int(len(self.data.sigma_cia_dict['xsecarr']))
        self.ctypes_cia_idx = (C.c_int*len(self.atmosphere.cia_idx))()
        self.ctypes_cia_idx[:] = self.atmosphere.cia_idx
        self.ctypes_cia_nidx = C.c_int(len(self.atmosphere.cia_idx))
        self.ctypes_sigma_cia_temp = cast2cpp(self.data.sigma_cia_dict['t'])
        self.ctypes_sigma_cia_ntemp = C.c_int(len(self.data.sigma_cia_dict['t']))
        self.ctypes_sigma_rayleigh_array = cast2cpp(self.atmosphere.sigma_rayleigh_array_flat)
        self.ctypes_clouds = C.c_int(1) if self.params.atm_clouds else C.c_int(0)
        self.ctypes_sigma_clouds_array = cast2cpp(self.atmosphere.sigma_clouds_array_flat)
        self.ctypes_star_radius = C.c_double(self.params.star_radius)

    def ctypes_update_vars(self):

        # load variables that will change during fitting

        self.ctypes_z = cast2cpp(self.atmosphere.altitude_profile)
        # dz = np.diff(self.atmosphere.altitude_profile)
        # self.ctypes_dz = cast2cpp(np.append(dz, dz[-1]))
        self.ctypes_active_mixratio_profile = cast2cpp(self.atmosphere.active_mixratio_profile)
        self.ctypes_inactive_mixratio_profile = cast2cpp(self.atmosphere.inactive_mixratio_profile)
        self.ctypes_temperature_profile = cast2cpp(self.atmosphere.temperature_profile)
        self.ctypes_planet_radius = C.c_double(self.atmosphere.planet_radius)
        self.ctypes_density_profile = cast2cpp(self.atmosphere.density_profile)
        self.ctypes_clouds_density_profile = cast2cpp(self.atmosphere.clouds_density_profile)

    #@profile
    def ctypes_pathintegral(self):

        self.ctypes_update_vars()
        #setting up output array
        absorption = zeros((self.atmosphere.int_nwngrid),dtype=float64)
        
        #retrieving function from cpp library
        ctypes_integral = self.cpathlib.path_integral

        #running c++ path integral
        ctypes_integral(self.ctypes_nwngrid,
                       self.ctypes_nlayers,
                       self.ctypes_nactive,
                       self.ctypes_ninactive,
                       self.ctypes_sigma_array,
                       self.ctypes_sigma_temp,
                       self.ctypes_sigma_ntemp,
                       self.ctypes_sigma_rayleigh_array,
                       self.ctypes_cia_npairs,
                       self.ctypes_cia_idx,
                       self.ctypes_cia_nidx,
                       self.ctypes_sigma_cia_array,
                       self.ctypes_sigma_cia_temp,
                       self.ctypes_sigma_cia_ntemp,
                       self.ctypes_clouds,
                       self.ctypes_sigma_clouds_array,
                       self.ctypes_clouds_density_profile,
                       self.ctypes_density_profile,
                       self.ctypes_z,
                       self.ctypes_active_mixratio_profile,
                       self.ctypes_inactive_mixratio_profile,
                       self.ctypes_temperature_profile,
                       self.ctypes_planet_radius,
                       self.ctypes_star_radius,
                       C.c_void_p(absorption.ctypes.data))

        out = np.zeros((len(absorption)))
        out[:] = absorption

        del(absorption)

        return out