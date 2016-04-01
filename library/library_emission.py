import scipy,itertools,sys,os,time,logging
import scipy.constants as con
import numpy as np
import ctypes as C
from scipy.stats.mstats_basic import tmean

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

def fit_brightness_temp(wave,flux):
    '''
    function fitting a black body to given flux-wavelength 
    value. 
    Input: Flux at wavelength (microns)
    Output: brightness temperature.
    '''
    
    def chi(temp,flux,wave):
        model = black_body(wave, temp)
        res = flux - model
        return res**2
    
    tempfit = scipy.optimize.fmin(chi,1000,args=(flux,wave),maxiter=1000,disp=0)
    return tempfit

def black_body_to_temp(wave,flux):
    '''
    Does the opposite to black_body. Converts flux to temperature. 
    Input: wavelength grid, flux grid 
    Output: tempeatrure grid
    @todo check if conversion is correct. There is probably a /m^2/micron offset 
    '''
    h = 6.62606957e-34
    c = 299792458
    k = 1.3806488e-23
    pi= 3.14159265359
    
    wave *= 1e-6
    
    logpart = np.log(((2.0*h*c**2)/(flux*wave**5))+1.0)
    T = (h*c)/(wave*k) * 1.0/ logpart
    return T

def iterate_TP_profile(TP_params, TP_params_std, TP_bounds, TP_function,iterate=True):
    '''
    function iterating through all lower and upper bounds of parameters
    to determine which combination gives the lowest/highest attainable
    TP profile. Returns mean TP profile with errorbars on each pressure level
    '''
    Tmean = TP_function(TP_params)

    bounds = [] #list of lower and upper parameter bounds
    lowpar = []
    highpar= []
    for i in xrange(len(TP_params)):
        low = TP_params[i]-TP_params_std[i]
        high = TP_params[i]+TP_params_std[i]
        lowpar.append(low)
        highpar.append(high)
        if low < TP_bounds[i][0]:
            low = TP_bounds[i][0]+1e-10
        if high > TP_bounds[i][1]:
            high = TP_bounds[i][1]-1e-10

        bounds.append((low,high))


    if iterate:
        iterlist = list(itertools.product(*bounds))
        iter_num = np.shape(iterlist)[0] #number of possible combinations

        T_iter   = np.zeros((len(Tmean), iter_num))
        T_minmax = np.zeros((len(Tmean), 2))

        for i in range(iter_num):
            T_iter[:,i]  = TP_function(iterlist[i])

        Tmean = np.mean(T_iter,1)
        T_minmax[:,0] = np.min(T_iter,1)
        T_minmax[:,1] = np.max(T_iter,1)
        T_sigma = (T_minmax[:,1] - T_minmax[:,0])/2.0
#         T_sigma = np.std(T_iter,1)

    else:
        Tmin = TP_function(lowpar)
        Tmax = TP_function(highpar)
        T_sigma = Tmax-Tmin
        T_sigma /= 2.0

    return Tmean, T_sigma


def generate_tp_covariance(outob):
    '''
    Function generating TP_profile covariance matrix from previous best fit. 
    This can be used by _TP_rodgers200 or _TP_hybrid TP profiles in a second stage fit.
    '''

    # todo needs to be adapted to new output class

    #translating fitting parameters to mean temperature and lower/upper bounds
    fit_TPparam_bounds = outob.fitting.fit_bounds[outob.fitting.fit_X_nparams:]
    if outob.NEST:
        T_mean, T_sigma = iterate_TP_profile(outob.NEST_TP_params_values[0], outob.NEST_TP_params_std[0],fit_TPparam_bounds,
                                              outob.fitting.forwardmodel.atmosphere.TP_profile)
    elif outob.MCMC:
        T_mean, T_sigma = iterate_TP_profile(outob.MCMC_TP_params_values[0], outob.MCMC_TP_params_std[0],fit_TPparam_bounds,
                                              outob.fitting.forwardmodel.atmosphere.TP_profile)
    elif outob.DOWN:
        FIT_std = np.zeros_like(outob.DOWN_TP_params_values)

        T_mean, T_sigma  = iterate_TP_profile(outob.DOWN_TP_params_values, FIT_std,fit_TPparam_bounds,
                                               outob.fitting.forwardmodel.atmosphere.TP_profile)
    else:
        logging.error('Cannot compute TP-covariance. No Stage 0 fit (NS/MCMC/MLE) can be found.')
        exit()
    
    #getting temperature error
    nlayers = outob.fitting.forwardmodel.atmosphere.nlayers
    
    #setting up arrays
    Ent_arr = np.zeros((nlayers,nlayers))
    Sig_arr = np.zeros((nlayers,nlayers))

    #populating arrays
    for i in range(nlayers):
        Ent_arr[i,:] = np.abs((T_mean[i])-(T_mean[:]))
        Sig_arr[i,:] = np.abs(T_sigma[i]+T_sigma[:])

    Diff_arr = np.sqrt(Ent_arr**2+Sig_arr**2)

    Diff_norm = ((Diff_arr-np.min(Diff_arr))/np.max(Diff_arr-np.min(Diff_arr)))
    Cov_array = 1.0 - Diff_norm
    
    return Cov_array