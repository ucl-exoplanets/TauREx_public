from numpy import *
# from pylab import *
import scipy
from scipy.special import wofz
# from library.library_occultquad import *
# from scipy.optimize import fmin
import sys
import os
# import string
# import math
import time
import ctypes as C

def interp_value(value, low_bound, up_bound, sig1, sig2):
    #interpolating a single value between low_bound and up_bound to new 
    #bounds sig1 and sig2. Probably quicker than np.interp function
    
    factor = (value - low_bound)/(up_bound - low_bound)
    
    new_value = sig1 + ((sig2-sig1)*factor)
    
    return new_value



