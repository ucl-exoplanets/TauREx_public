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
        if self.params.gen_trans_cpp:
            self.cpathlib = C.CDLL('./library/cpath_pathintegral.so', mode=C.RTLD_GLOBAL)

        # preload variables in memory for cpath
        self.cpath_load_vars()

        # set forward model function
        self.model = self.cpath_pathintegral

    # class methods

    def get_cloud_sig(self):
        # calculating could cross sections
        a = (128.0* pi**5 * self.atmosphere.clouds_a**6)
        b = (3.0 * (10000./self.data.obs_wngrid)**4)
        c = ((self.atmosphere.clouds_m**2 -1.0)/(self.atmosphere.clouds_m**2 + 2.0))**2
        return a / b * c

    def cpath_load_vars(self):

        # load variables that won't change during fitting
        self.cpath_nwngrid = C.c_int(self.data.int_nwngrid)
        self.cpath_nlayers = C.c_int(self.atmosphere.nlayers)
        self.cpath_nactive = C.c_int(self.atmosphere.nactivegases)
        self.cpath_ninactive = C.c_int(self.atmosphere.ninactivegases)
        self.cpath_sigma_array = cast2cpp(self.atmosphere.sigma_array_flat)
        self.cpath_sigma_temp = cast2cpp(self.data.sigma_dict['t'])
        self.cpath_sigma_ntemp = C.c_int(len(self.data.sigma_dict['t']))
        self.cpath_sigma_rayleigh_array = cast2cpp(self.atmosphere.sigma_rayleigh_array_flat)
        self.cpath_cia_npairs = C.c_int(len(self.data.sigma_cia_dict['xsecarr']))
        self.cpath_cia_idx = (C.c_int*len(self.atmosphere.cia_idx))()
        self.cpath_cia_idx[:] = self.atmosphere.cia_idx
        self.cpath_cia_nidx = C.c_int(len(self.atmosphere.cia_idx))
        self.cpath_sigma_cia_array = cast2cpp(self.atmosphere.sigma_cia_array_flat)
        self.cpath_sigma_cia_temp = cast2cpp(self.data.sigma_cia_dict['t'])
        self.cpath_sigma_cia_ntemp = C.c_int(len(self.data.sigma_cia_dict['t']))
        self.cpath_star_radius = C.c_double(self.params.star_radius)

    def cpath_update_vars(self):

        # load variables that will change during fitting
        self.cpath_z = cast2cpp(self.atmosphere.altitude_profile)
        dz = np.diff(self.atmosphere.altitude_profile)
        self.cpath_dz = cast2cpp(np.append(dz, dz[-1]))
        self.cpath_density_profile = cast2cpp(self.atmosphere.density_profile)
        self.cpath_active_mixratio_profile = cast2cpp(self.atmosphere.active_mixratio_profile)
        self.cpath_inactive_mixratio_profile = cast2cpp(self.atmosphere.inactive_mixratio_profile)
        self.cpath_temperature_profile = cast2cpp(self.atmosphere.temperature_profile)
        self.cpath_planet_radius = C.c_double(self.atmosphere.planet_radius)

    #@profile
    def cpath_pathintegral(self):

        self.cpath_update_vars()
    
        #setting up output array
        absorption = zeros((self.data.int_nwngrid),dtype=float64)
        
        #retrieving function from cpp library
        cpath_integral = self.cpathlib.path_integral

        #running c++ path integral
        cpath_integral(self.cpath_nwngrid,
                       self.cpath_nlayers,
                       self.cpath_nactive,
                       self.cpath_ninactive,
                       self.cpath_sigma_array,
                       self.cpath_sigma_temp,
                       self.cpath_sigma_ntemp,
                       self.cpath_sigma_rayleigh_array,
                       self.cpath_cia_npairs,
                       self.cpath_cia_idx,
                       self.cpath_cia_nidx,
                       self.cpath_sigma_cia_array,
                       self.cpath_sigma_cia_temp,
                       self.cpath_sigma_cia_ntemp,
                       self.cpath_z,
                       self.cpath_dz,
                       self.cpath_density_profile,
                       self.cpath_active_mixratio_profile,
                       self.cpath_inactive_mixratio_profile,
                       self.cpath_temperature_profile,
                       self.cpath_planet_radius,
                       self.cpath_star_radius,
                       C.c_void_p(absorption.ctypes.data))

        out = np.zeros((len(absorption)))
        out[:] = absorption

        del(absorption)

        return out