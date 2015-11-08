from numpy import *
from numpy.ctypeslib import as_ctypes
# from library.library_preselector import *
import scipy, sys, glob, os, time,string
from scipy.special import wofz
import ctypes as C
# from pylab import *
# from library.library_occultquad import *
# from scipy.optimize import fmin
# import string
import numpy as np
import math

# def linear_abs_temp_interpolator(TEMPLIST, ):

def house_keeping(params,options):
    #does some housekeeping for the final fitting results

    #copies used parameter file to ./Output
    if params.clean_save_used_params:
        subprocess.call('cp '+options.param_filename+' Output/',shell=True)

    subprocess.call('python '+params.clean_script,shell=True)


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

def convert2microns(path, upcut=30):
#Function converting ExoMol cross section files in dir:PATH from wavenumbers to microns and sorting
#with ascending wavelength
#input: .txt files in wavenumbers
#output: .abs files in microns

    for path, subdirs, files in os.walk(path):
        # list files in all subdirectories recursively
        for name in files:
            filename = os.path.join(path, name)
            if filename[-3:] ==  'txt':
                if os.path.isfile(filename[:-3] + 'abs'):
                    # skip file if it already exists
                    pass
                else:
                    print 'Converting: %s' % filename
                    tmp = loadtxt(filename, dtype=float32)[1:,:]
                    tmp[:,0] = 10000.0/tmp[:,0]
                    idx = argsort(tmp[:,0],axis=-1)
                    tmp2 = tmp[idx,:][where(tmp[idx,0] < upcut)]
                    tmp2 = tmp2.astype(float32, copy=False)
                    savetxt(filename[:-3]+'abs', tmp2, fmt="%.6e %.8e")

def find_nearest(arr, value):
    # find nearest value in array
    arr = array(arr)
    idx = (abs(arr-value)).argmin()
    return [arr[idx], idx]


def find_absfiles(path, molname, ext='sigma'):
    # finding all absorption crosssection files in path for molecule molname (path/molname/*.ext)
    # filename convection follows that of zero pressure ExoMol cross sections
    # return: array of absfilenames and array of corresponding temperatures

    globlist = glob.glob(os.path.join(path, molname, '*.%s' % ext))

    absfilelist = []
    templist = []

    for file in globlist:
        fname = string.rsplit(file,'/',1)[1] #splitting the name
        splitname = string.split(fname,'_',3)
        absfilelist.append(fname)
        templist.append(float(splitname[2][:-1]))

    templist = asarray(templist)
    absfilelist = asarray(absfilelist)

    sortidx = argsort(templist)

    return absfilelist[sortidx], templist[sortidx]

def find_absfiles_pressure(path, molname):
    # finding all absorption crosssection files in path for molecule molname (path/molname/molanme_*.abs)
    # filename convection follows that of first generation of pressure broadened ExoMol cross sections
    # return: array of absfilenames and array of corresponding temperatures and pressures

    absfilelist = {}
    templist = []
    presslist = []

    for file in glob.glob(os.path.join(path, molname, '*.sig')):
        fname = os.path.basename(file)
        splitname = string.split(fname,'_',5)
        temp = float(splitname[3])
        pres = float(splitname[4])
        templist.append(temp)
        presslist.append(pres)
        if not temp in absfilelist:
            absfilelist[temp] = {}
        absfilelist[temp][pres] = fname

    return absfilelist, list(set(templist)), list(set(presslist))



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
#         ARRAY = ARRAY.astype(float64)
  
    ARRdim = len(shape(ARRAY)) #getting number of dimensions
    if ARRdim == 1:
        s1 = len(ARRAY)
        return ARRAY.ctypes.data_as(C.POINTER(C.c_double)) #creating 1D pointer array
    elif ARRdim == 2:
        [s1,s2] = shape(ARRAY)
        dbptr = C.POINTER(C.c_double)
        PARR = (dbptr*s1)(*[row.ctypes.data_as(dbptr) for row in ARRAY]) #creating 2D pointer array
        return PARR
    elif ARRdim == 3:
        PARR = ARRAY.ctypes.data_as(C.POINTER(C.c_double))
        [s1,s2,s3] = shape(ARRAY)
        dbptr = C.POINTER(C.c_double)
        PARR = (dbptr*s1)(dbptr*s2)(*[row.ctypes.data_as(dbptr) for row in ARRAY]) #creating 2D pointer array
#         PARR = as_ctypes(ARRAY)
        
        return PARR
   

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

def round_base(x, base=1):
    return int(base * round(float(x)/base))

def binspectrum(spectrum_in, resolution):
    wavegrid, dlamb_grid = get_specgrid(R=resolution,lambda_min=np.min(spectrum_in[:,0]),lambda_max=np.max(spectrum_in[:,0]))
    spec_bin_grid, spec_bin_grid_idx = get_specbingrid(wavegrid, spectrum_in[:,0])
    spectrum_binned = [spectrum_in[:,1][spec_bin_grid_idx == i].mean() for i in xrange(1,len(spec_bin_grid))]
    return transpose(vstack((spec_bin_grid[:-1], spectrum_binned)))

def get_specgrid( R=5000, lambda_min=0.1, lambda_max=20.0):
    #generating wavelength grid with uniform binning in log(lambda)
    #lambda min and max are in microns, R is the spectral resolution
    #R = lambda/delta_lamdba
    specgrid = []
    delta_lambda =[]
    specgrid.append(lambda_min)
    run = True
    i=1
    while run:
        dlam= specgrid[i-1]/R
        specgrid.append(specgrid[i-1]+dlam)
        delta_lambda.append(dlam)

        if specgrid[i] >= lambda_max:
            run=False
        i+=1
    return np.asarray(specgrid),np.asarray(delta_lambda)

def get_specbingrid(wavegrid, specgrid, binwidths=None):
    #function calculating the bin boundaries for the data
    #this is used to bin the internal spectrum to the data in fitting module
    if not isinstance(binwidths, (np.ndarray, np.generic)):
        bingrid =[]
        bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
        for i in range(len(wavegrid)-1):
            bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
        bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
        bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
    else:
        # this bingrid is actually useless, as it doesn't allow for gaps in the data
        bingrid = []
        for i in range(len(wavegrid)):
            bingrid.append(wavegrid[i]-binwidths[i]/2.)
            bingrid.append(wavegrid[i]+binwidths[i]/2.)

        # build bin grid index array (an index for each model datapoint
        bingrid_idx = np.empty(len(specgrid))
        bingrid_idx[:] = np.NaN
        for i in range(len(specgrid)):
            for j in range(len(wavegrid)):
                if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                    bingrid_idx[i] = j+1
                    break
    return bingrid, bingrid_idx

def plot_bin(spectrum, R, ycol=1, yadg=0., **kwargs):
    sptmp = np.zeros((len(spectrum[:,0]),2))
    sptmp[:,0] = spectrum[:,0]
    sptmp[:,1] = spectrum[:,ycol]
    spectrum_bin = binspectrum(sptmp, R)
    plt.plot(spectrum_bin[:,0], spectrum_bin[:,1]+yadg, **kwargs)


