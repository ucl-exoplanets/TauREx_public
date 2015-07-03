import sys
cimport cython

# sys.path.append('./library')
import library_emission as em
import library_general as gen


import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

@cython.boundscheck(False) # turn of bounds-checking for entire function
cdef black_body(np.ndarray[DTYPE_t,ndim=1] lamb, int nlambda, np.float_t temp):
    #small function calculating plank black body
    #input: microns, kelvin
    #output: W/m^2/micron
    
    cdef DTYPE_t h, c, k, pi
    cdef Py_ssize_t l
    
    cdef np.ndarray[DTYPE_t,ndim=1] BB = np.zeros([nlambda],dtype=DTYPE)
    
    h = 6.62606957e-34
    c = 299792458
    k = 1.3806488e-23
    pi= 3.14159265359
    

    BB = 1e-6 * (pi* (4.0*h*c**2)/(lamb*1e-6)**5) * (1.0/(np.exp((h * c) / (lamb *1e-6 * k * temp)) -1))
    return BB

@cython.boundscheck(False) # turn of bounds-checking for entire function
cdef np.ndarray[DTYPE_t,ndim=2] get_sigma_array(dict sigma_dict,np.float_t temperature):
    #getting sigma array from sigma_dic for given temperature
    return sigma_dict[gen.find_nearest(sigma_dict['tempgrid'],temperature)[0]]

@cython.boundscheck(False) # turn of bounds-checking for entire function
def path_integral(np.ndarray[DTYPE_t,ndim=2] X, np.ndarray[DTYPE_t,ndim=1] rho, np.ndarray[DTYPE_t,ndim=1] temperature, dict sigma_dict,np.ndarray[DTYPE_t,ndim=1] BB_star, 
                  np.ndarray[DTYPE_t,ndim=1] specgrid, np.ndarray[DTYPE_t,ndim=1] dzarray, int nlambda, int nlayers, int n_gas):
    
    #array definitions 
    cdef np.ndarray[DTYPE_t,ndim=1] I_total   = np.zeros([nlambda],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] tau       = np.zeros([nlayers,nlambda],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] dtau      = np.zeros([nlayers,nlambda],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=2] tau_total = np.zeros([nlayers,nlambda],dtype=DTYPE)
    
    cdef np.ndarray[DTYPE_t,ndim=2] sigma_array  = np.zeros([n_gas,nlambda],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] BB_surf      = np.zeros([nlambda],dtype=DTYPE)
    cdef np.ndarray[DTYPE_t,ndim=1] BB_layer     = np.zeros([nlambda],dtype=DTYPE)
    
    cdef Py_ssize_t k, i, j, l
    
    #surface layer      
    BB_surf = black_body(specgrid,nlambda,temperature[0])  
#     BB_surf = em.black_body(specgrid,temperature[0])
    sigma_array = get_sigma_array(sigma_dict,temperature[0])

    for k in xrange(nlayers):
        if temperature[k] != temperature[k-1]:
            sigma_array = get_sigma_array(sigma_dict,temperature[k])
                
        for i in xrange(n_gas):
            for l in xrange(nlambda):
                tau[0,l] += (sigma_array[i,l] * X[i,k] * rho[k] * dzarray[k])
                
    I_total += BB_surf*(np.exp(-1.0*tau[0,:]))
        
    #other layers
    BB_layer = BB_surf
    
    sigma_array = get_sigma_array(sigma_dict,temperature[0])
    for j in xrange(1,nlayers):   
        if temperature[j] != temperature[j-1]: 
            BB_layer = black_body(specgrid,nlambda,temperature[j]) 
#             BB_layer = em.black_body(specgrid,temperature[j])
                        
        for k in xrange(j,nlayers):
            if temperature[k] != temperature[k-1]:
                sigma_array = get_sigma_array(sigma_dict,temperature[k])  
                
            if j is k:
                for i in xrange(n_gas):
                    for l in xrange(nlambda):
                        tau[j,l] += (sigma_array[i,l] * X[i,k] * rho[k] * dzarray[k])
                        dtau[j,l] += (sigma_array[i,l] * X[i,j] * rho[j] * dzarray[j])
            else:
                for i in xrange(n_gas):
                    for l in xrange(nlambda):
                        tau[j,l] += (sigma_array[i,l] * X[i,k] * rho[k] * dzarray[k])         
          
        
        tau_total[j,:] = BB_layer*(np.exp(-1.0*tau[j,:])) * (dtau[j,:])
        I_total += BB_layer*(np.exp(-1.0*tau[j,:])) * (dtau[j,:])

    return I_total, tau_total