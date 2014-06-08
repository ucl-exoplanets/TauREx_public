from numpy import *
# from library.library_preselector import *
import scipy, sys, glob, os, time,string
from scipy.special import wofz
import ctypes as C
# from pylab import *
# from library.library_occultquad import *
# from scipy.optimize import fmin
# import string
# import math

# def linear_abs_temp_interpolator(TEMPLIST, ):


def find_nearest(arr, value):
    # find nearest value in array
    arr = array(arr)
    idx = (abs(arr-value)).argmin()
    return [arr[idx], idx]


def find_absfiles(PATH, MOLNAME):
    #finding all absorption crosssection files in PATH for molecule MOLNAME
    #return: array of absfilenames and array of corresponding temperatures
    
    globlist = glob.glob(PATH+'*.abs')
    
    absfilelist = []
    templist = []

    for FILE in globlist:
        fname = string.rsplit(FILE,'/',1)[1] #splitting the name
        splitname = string.split(fname,'_',3)

        if splitname[0] == MOLNAME:
            absfilelist.append(fname)
            templist.append(float(splitname[2][:-1]))
    
    templist = asarray(templist)
    absfilelist = asarray(absfilelist)
    
    sortidx = argsort(templist)
    
#     print absfilelist[sortidx], templist[sortidx]
    
    return absfilelist[sortidx], templist[sortidx]



def find_single_absfile(PATH,MOLECULELIST,TEMPERATURE):
    #determining correct abs files to be read in
    #reading available cross section lists in PATH
    
    #code can be improved
    
    globlist = glob.glob(PATH+'*.abs')

    #determining the right abs file for correct Tplanet
    #this needs to be changed when we want several temperatures

    absfilelist = []
    mollist = []
    
    #first pass looking for available molecules
    for FILE in globlist:
        fname = string.rsplit(FILE,'/',1)[1] #splitting the name
        splitname = string.split(fname,'_',3)
        
        if splitname[0] not in mollist:
            mollist.append(splitname[0])
            
            
    #second pass looking for available temperatures
    for molecule in MOLECULELIST:
        if molecule not in mollist:
            raise ValueError('Absorption cross-section file not found for: '+molecule)
        else:
            temps = []
            for FILE in globlist:
                fname = string.rsplit(FILE,'/',1)[1] #splitting the name
                splitname = string.split(fname,'_',3)
                temps.append(float(splitname[2][:-1])) #getting temperature from file name
                
            #looking for closest temperature match
            next_temp = find_nearest(temps,TEMPERATURE)[0]
        
            #selecting correct cross section file
            for FILE in globlist:
                fname = string.rsplit(FILE,'/',1)[1] #splitting the name
                splitname = string.split(fname,'_',3)
                
                if splitname[0] == molecule:
                    if float(splitname[2][:-1]) == next_temp:
                        absfilelist.append(fname)  
    
    return absfilelist




def cast2cpp(ARRAY):
    #cast numpy array to C array

    if ARRAY.dtype != float64:
        print 'WARNING: array not cast in float64'
        print 'current format: ', ARRAY.dtype

    ARRdim = len(shape(ARRAY)) #getting number of dimensions
    if ARRdim == 1:
        s1 = len(ARRAY)
        return ARRAY.ctypes.data_as(C.POINTER(C.c_double)), C.c_int(s1) #creating 1D pointer array
    elif ARRdim == 2:
        [s1,s2] = shape(ARRAY)
        dbptr = C.POINTER(C.c_double)
        PARR = (dbptr*s1)(*[row.ctypes.data_as(dbptr) for row in ARRAY]) #creating 2D pointer array
        return PARR, C.c_int(s1), C.c_int(s2)
    
def redirect_stderr_stdout(stderr=sys.stderr, stdout=sys.stdout):
    #function decorator to re-direct standard output
    #use as: @redirect_stderr_stdout(some_logging_stream, the_console)
    #to print to file: @redirect_stderr_stdout(stdout=open(FILENAME,'wb'))
    def wrap(f):
        def newf(*args, **kwargs):
            old_stderr, old_stdout = sys.stderr, sys.stdout
            sys.stderr = stderr
            sys.stdout = stdout
            try:
                return f(*args, **kwargs)
            finally:
                sys.stderr, sys.stdout = old_stderr, old_stdout

        return newf
    return wrap

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "="*block + " "*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()   
    
    
    
def voigt(x, y):
    # The Voigt function is also the real part of 
    # w(z) = exp(-z^2) erfc(iz), the complex probability function,
    # which is also known as the Faddeeva function. Scipy has 
    # implemented this function under the name wofz()
    z = x + 1j*y
    I = wofz(z).real
    return I

def Voigt(nu, alphaD, alphaL, nu_0, A, a=0, b=0):
    # The Voigt line shape in terms of its physical parameters
    f = sqrt(log(2))
    x = (nu-nu_0)/alphaD * f
    y = alphaL/alphaD * f
    backg = a + b*nu 
    V = A*f/(alphaD*sqrt(pi)) * voigt(x, y) + backg
    return V

def funcVoigt(p, x):
    # Compose the Voigt line-shape
    alphaD, alphaL, nu_0, I, a, b = p
    return Voigt(x, alphaD, alphaL, nu_0, I, a, b)

def funcGauss(p, x):
    # Gaussian function
    mu = p[0]
    sigma = p[1]
    return 1.0/(sigma * sqrt(2.0 * pi)) * exp(-(x-mu)**2/(2.0*sigma**2))


        