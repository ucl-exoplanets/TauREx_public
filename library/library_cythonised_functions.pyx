#library containing cythonised functions for added speed
#the library can be compiled from taurex root dir with: 
# python cython_setup.py build_ext --inplace  
# or 
# python c_setup.py build_ext --inplace
#
# For functions to be callable, this library does  need to be compiled 
#
###########################

cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
def get_sigma_array_interp(
    np.ndarray[DTYPE_t,ndim=1] temp_list,
    np.ndarray[DTYPE_t,ndim=1] press_list,
    pressure_val,
    np.ndarray[DTYPE_t,ndim=3] sigma_in_cut,
    np.ndarray[DTYPE_t,ndim=1] int_wngrid,
    int_wngrid_idxmin,
    int_wngrid_idxmax):

    cdef np.ndarray[DTYPE_t,ndim=2] sigma_array_partial = np.zeros((len(temp_list), len(int_wngrid)),dtype=DTYPE)
    cdef Py_ssize_t i # ??

    for temperature_idx in range(len(temp_list)):
        for wno_idx in range(len(int_wngrid)):
            sigma_array_partial[temperature_idx, wno_idx] = \
                np.interp(pressure_val, press_list, sigma_in_cut[:,temperature_idx,wno_idx])

    return sigma_array_partial

@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
def runtime_bin_spectrum(np.ndarray[DTYPE_t,ndim=1] model,np.ndarray[np.int_t,ndim=1] spec_bin_grid_idx,int n_spec_bin_grid):

    cdef np.ndarray[DTYPE_t,ndim=1] model_binned = np.zeros((n_spec_bin_grid),dtype=DTYPE)
    cdef Py_ssize_t i

    for i in range(1, n_spec_bin_grid+1):
        model_binned[i] = np.mean(model[spec_bin_grid_idx == i]) #@todo this is not faster than the python version. need to come up with something more clever...

    return model_binned
