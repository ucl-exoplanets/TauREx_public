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
def runtime_bin_spectrum(np.ndarray[DTYPE_t,ndim=1] model,np.ndarray[np.int_t,ndim=1] spec_bin_grid_idx,int n_spec_bin_grid):
    
    cdef np.ndarray[DTYPE_t,ndim=1] model_binned = np.zeros((n_spec_bin_grid),dtype=DTYPE)
    cdef Py_ssize_t i
    
    for i in xrange(1, n_spec_bin_grid+1):
        model_binned[i] = np.mean(model[spec_bin_grid_idx == i]) #@todo this is not faster than the python version. need to come up with something more clever...
        
    return model_binned