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

    def __init__(self, atmosphere, stage=0, data=None, params=None):

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

        self.stage = 0

        if self.params.gen_ace:
            # loading Fortran code for chemically consistent model
            self.ace_lib = C.CDLL('./library/ACE/ACE.so', mode=C.RTLD_GLOBAL)

        if self.params.in_opacity_method in ['xsec_highres', 'xsec_lowres']: # using cross sections

            if self.atmosphere.nthreads > 1:
                # load openmp version. Remember to set OMP_NUM_THREADS to number of threads
                self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_emission_parallel.so', mode=C.RTLD_GLOBAL)
            else:
                self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_emission.so', mode=C.RTLD_GLOBAL)

            # set forward model function
            self.model = self.ctypes_pathintegral

            #retrieving function from cpp library
            self.pathintegral_lib.path_integral.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
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
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 C.c_double,
                 C.c_double,
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 C.c_void_p]



        elif self.params.in_opacity_method in ['ktab', 'ktable', 'ktables']: # using k tables

            self.pathintegral_lib = C.CDLL('./library/ctypes_pathintegral_emission_ktab.so', mode=C.RTLD_GLOBAL)
            self.model = self.ctypes_pathintegral_ktab

            #retrieving function from cpp library
            self.pathintegral_lib.path_integral.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 C.c_int,
                 C.c_int,
                 C.c_int,
                 C.c_int,
                 C.c_int,
                 C.c_int,
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
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 C.c_double,
                 C.c_double,
                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
                 C.c_void_p]

    def ctypes_pathintegral(self, return_tau=False, mixratio_mask=False):
        
        if self.params.gen_ace:
            self.set_ACE(mixratio_mask)

        #setting up output array
        FpFs = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_wngrid,
                                             self.atmosphere.int_nwngrid,
                                             self.atmosphere.nlayers,
                                             self.atmosphere.nactivegases,
                                             self.atmosphere.ninactivegases,
                                             self.params.atm_rayleigh,
                                             self.params.atm_cia,
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
                                             self.atmosphere.density_profile,
                                             self.atmosphere.altitude_profile,
                                             self.atmosphere.active_mixratio_profile.flatten(),
                                             self.atmosphere.inactive_mixratio_profile.flatten(),
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
            for i in range(self.atmosphere.int_nwngrid):
                for j in range(self.atmosphere.nlayers-1):
                    tauout[i,j] = tau[i + j*self.atmosphere.int_nwngrid]

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

    def ctypes_pathintegral_ktab(self, return_tau=False, mixratio_mask=False):


        if self.params.gen_ace:

            self.set_ACE(mixratio_mask)

        #setting up output array
        FpFs = zeros((self.atmosphere.int_nwngrid), dtype=np.float64, order='C')
        tau = zeros((self.atmosphere.int_nwngrid*self.atmosphere.nlayers), dtype=np.float64, order='C')

        #running c++ path integral
        self.pathintegral_lib.path_integral(self.atmosphere.int_wngrid,
                                             self.atmosphere.int_nwngrid,
                                             self.atmosphere.nlayers,
                                             self.atmosphere.nactivegases,
                                             self.atmosphere.ninactivegases,
                                             self.params.atm_rayleigh,
                                             self.params.atm_cia,
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
                                             self.atmosphere.density_profile,
                                             self.atmosphere.altitude_profile,
                                             self.atmosphere.active_mixratio_profile.flatten(),
                                             self.atmosphere.inactive_mixratio_profile.flatten(),
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
            for i in range(self.atmosphere.int_nwngrid):
                for j in range(self.atmosphere.nlayers-1):
                    tauout[i,j] = tau[i + j*self.atmosphere.int_nwngrid]

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


    def set_ACE(self, mixratio_mask):

        # todo should really move this to atmosphere.py [see also transmission]

        # chemically consistent model
        vector = C.c_double*self.atmosphere.nlayers
        a_apt = vector()
        p_apt = vector()
        t_apt = vector()
        for i in range(self.atmosphere.nlayers):
           a_apt[i] = self.atmosphere.altitude_profile[i]/1000.
           p_apt[i] = self.atmosphere.pressure_profile[i]/1.e5
           t_apt[i] = self.atmosphere.temperature_profile[i]

        # y_out has shape (nlayers, 105). 105 is the total number of molecules computed
        y_out = ((C.c_double * 105) * self.atmosphere.nlayers)()

        self.ace_lib.ACE(C.byref(C.c_int(self.atmosphere.nlayers)),
                         C.byref(a_apt),
                         C.byref(p_apt),
                         C.byref(t_apt),
                         C.byref(C.c_double(self.atmosphere.He_abund_dex)),
                         C.byref(C.c_double(self.atmosphere.C_abund_dex)),
                         C.byref(C.c_double(self.atmosphere.O_abund_dex)),
                         C.byref(C.c_double(self.atmosphere.N_abund_dex)),
                         C.byref(y_out))
        ace_profiles = np.asarray(y_out)


        for mol_idx, mol_val in enumerate(self.params.atm_active_gases):
            self.atmosphere.active_mixratio_profile[mol_idx, :] = ace_profiles[:, self.data.ace_active_gases_idx[mol_idx]]

        for mol_idx, mol_val in enumerate(self.params.atm_inactive_gases):
            self.atmosphere.inactive_mixratio_profile[mol_idx, :] = ace_profiles[:, self.data.ace_inactive_gases_idx[mol_idx]]

        if isinstance(mixratio_mask, (np.ndarray, np.generic)):
            self.atmosphere.active_mixratio_profile[mixratio_mask, :] = 0

        del(y_out)
        del(a_apt)
        del(p_apt)
        del(t_apt)

        # couple mu to composition
        self.atmosphere.planet_mu = self.atmosphere.get_coupled_planet_mu()

        # update atmospheric params
        self.atmosphere.set_altitude_gravity_scaleheight_profile()


