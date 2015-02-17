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


def generate_tp_covariance(fitob):
    '''
    Function generating TP_profile covariance matrix from previous best fit. 
    This can be used by _TP_rodgers200 or _TP_hybrid TP profiles in a second stage fit.
    '''
    
    #translating fitting parameters to mean temperature and lower/upper bounds
    T_mean, T_minmax, P = emlib.iterate_TP_profile(fitob.NEST_FIT_mean[0], fitob.NEST_FIT_std[0], 
                                             fitob.atmosphere.fit_index,fitob.atmosphere.TP_profile)
    
    #getting temperature error
    T_sigma = T_minmax[:,1] - T_minmax[:,0]
    nlayers = fitob.atmosphere.nlayers
    
    #setting up arrays
    Ent_arr = np.zeros((nlayers,nlayers))
    Sig_arr = np.zeros((nlayers,nlayers))
    
    #populating arrays
    for i in range(nlayers):
            Ent_arr[i,:] = np.abs(T_mean[i]-T_mean[:])
            Sig_arr[i,:] = np.sqrt(T_sigma[i]**2+T_sigma[:]**2)

    Diff_arr = Ent_arr-Sig_arr
    Diff_norm = ((Diff_arr-np.min(Diff_arr))/np.max(Diff_arr-np.min(Diff_arr)))
    Cov_array = 1.0 - Diff_norm
    
    return Cov_array