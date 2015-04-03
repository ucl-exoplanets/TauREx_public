

#loading libraries
import numpy, pylab, sys, os, optparse, time
from numpy import * #nummerical array library
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

#checking for multinest library
global multinest_import
try:
    with os.devnull as sys.stdout:
        import pymultinest
        multinest_import = True
except:
    multinest_import = False

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,atmosphere,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

## Apply TP profile
def movingaverage(values,window):
    weigths = np.repeat(1.0, window)/window
    smas = np.convolve(values, weigths, 'valid')
    return smas #smas2[::-1] # as a numpy array


AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  help='Parameter filename')
options, remainder = parser.parse_args()

#initialising parameters object
params = parameters(options.param_filename)

# set model resolution to 1000
params.gen_spec_res = 1000
params.gen_manual_waverange = True

#errors = [50, 100, 300]
#resolutions = [300, 100, 50, 10]

resolution = 1000

params.gen_manual_waverange = True
params.in_spectrum_file = False


#initialising data object
dataob = data(params)

#initialising TP profile object
atmosphereob = atmosphere(dataob)

mixing = [atmosphereob.absorbing_gases_X[idx] for idx, gasname in enumerate(atmosphereob.absorbing_gases_X)]


for j in range(len(atmosphereob.absorbing_gases) + 3):

    if j > 0:
        for i in range(len(atmosphereob.absorbing_gases_X)):
            atmosphereob.absorbing_gases_X[i] = 1e-20

    if j == 1: # rayleigh
        atmosphereob.params.in_include_cld = False
    elif j == 2: # clouds
        atmosphereob.params.in_include_cld = True
    elif j > 2:
        idx = j-3
        atmosphereob.params.in_include_cld = False
        atmosphereob.absorbing_gases_X[idx] = mixing[idx]

    atmosphereob.X = atmosphereob.set_mixing_ratios()

    logging.info('Creating custom TP profile')
    MAX_P = atmosphereob.P[0]
    MIN_P = atmosphereob.P[-1]
    smooth_window = 5 #smoothing window size as percent of total data
    Pnodes = [MAX_P, MAX_P, 0.65e4, MIN_P]
    Tnodes = [1000,1000, 600,600]
    TP = np.interp((np.log(atmosphereob.P[::-1])), np.log(Pnodes[::-1]), Tnodes[::-1])
    atmosphereob.T = TP[::-1]
    atmosphereob.update_atmosphere()

    # save TP profile to file
    # out = np.zeros((len(atmosphereob.T),2))
    # out[:,0] = atmosphereob.T
    # out[:,1] = atmosphereob.P
    # np.savetxt(os.path.join(params.out_path, 'TP_profile.dat'), out)

    #initialising transmission radiative transfer code object
    forwardmodelob = transmission(atmosphereob)

    # bin down internal model to given resolution
    wavegrid, dlamb_grid = dataob.get_specgrid(R=resolution,lambda_min=params.gen_wavemin,lambda_max=params.gen_wavemax)
    spec_bin_grid, spec_bin_grid_idx = dataob.get_specbingrid(wavegrid, dataob.specgrid)

    model = forwardmodelob.model()
    model_binned = [np.mean(model[spec_bin_grid_idx == i]) for i in xrange(1,len(spec_bin_grid))]

    out = np.zeros((len(wavegrid),3))
    out[:,0] = wavegrid
    out[:,1] = model_binned

    # remove nan values. Sometimes the last bin has a nan...
    out = out[~np.isnan(out).any(axis=1)]

    #initiating output object with fitted data from fitting class
    outputob = output(forwardmodel=forwardmodelob)

    if j == 0:
        outputob.save_spectrum_to_file(spectrum=out, saveas='all_spectrum.dat')
    elif j == 1: # rayleigh
        outputob.save_spectrum_to_file(spectrum=out, saveas='rayleigh_spectrum.dat')
    elif j == 2: # clouds
        outputob.save_spectrum_to_file(spectrum=out, saveas='clouds_spectrum.dat')
    else:
        outputob.save_spectrum_to_file(spectrum=out, saveas='%s_spectrum.dat' % atmosphereob.absorbing_gases[idx])
