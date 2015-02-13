import scipy,itertools,sys,os,time
import scipy.constants as con
import numpy as np
import ctypes as C


def black_body(lamb, temp):
    #small function calculating plank black body
    #input: microns, kelvin
    #output: W/m^2/micron
    
    h = 6.62606957e-34
    c = 299792458
    k = 1.3806488e-23
    pi= 3.14159265359
    
    exponent = np.exp((h * c) / (lamb *1e-6 * k * temp))
    BB = (pi* (2.0*h*c**2)/(lamb*1e-6)**5) * (1.0/(exponent -1))
    
#     exponent = np.exp((con.h * con.c) / (lamb *1e-6 * con.k * temp))
#     BB = (np.pi* (2.0*con.h*con.c**2)/(lamb*1e-6)**5) * (1.0/(exponent -1))
    
    return BB * 1e-6


def iterate_TP_profile(fit_params, fit_params_std, fit_idx, TP_function):
    '''
    function iterating through all lower and upper bounds of parameters
    to determine which combination gives the lowest/highest attainable 
    TP profile. Returns mean TP profile with errorbars on each pressure level
    ''' 
    
    TP_params     = fit_params[fit_idx[0]:]
    TP_params_std = fit_params_std[fit_idx[0]:]
    
    
    Tmean,P,X = TP_function(np.asarray(fit_params))
    
    bounds = [] #list of lower and upper parameter bounds 
    for i in xrange(len(TP_params)):
        bounds.append((TP_params[i]-TP_params_std[i],TP_params[i]+TP_params_std[i]))
    
    iterlist = list(itertools.product(*bounds))
    iter_num = np.shape(iterlist)[0] #number of possible combinations
    
 
    T_iter   = np.zeros((len(Tmean),iter_num))
    T_minmax = np.zeros((len(Tmean),2))
    
    for i in range(iter_num):
        iterlist2 = fit_params
        iterlist2[fit_idx[0]:] = list(iterlist[i])
        T,P, X = TP_function(iterlist2)
        T_iter[:,i] = T
    
    T_minmax[:,0] = np.min(T_iter,1)
    T_minmax[:,1] = np.max(T_iter,1)
    
    #Say hello, to the rug's topography
    #It holds quite a lot of interest with your face down on it. 
    #Say hello, to the shrinking in your head
    #You can't see it but you know it's there, so don't negelect it. 
    
    #I'm taking her home with me. All dressed in white.
    #She's got everything I need: pharmacy keys.
    #She's fallen hard for me. I can see it in her eyes
    #She acts just like a nurse... with all the other guys. 
    
    return Tmean, T_minmax, P