import numpy as np
# from pylab import *
import scipy
import scipy.constants as con
# from library.library_occultquad import *
# from scipy.optimize import fmin
import sys
import os
# import string
# import math
import time
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