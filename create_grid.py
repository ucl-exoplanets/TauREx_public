#! /usr/bin/python -W ignore

###########################################################
# create_grid
# create transmission spectra from grid of planets, spectral resolutions
# and signal to noise
#
# Requirements: -python libraries: pylab, numpy, ConfigParser
#             [these are the minimum requirements]
#
# Additional requirements: pymc
#
# Inputs: nothing
#
# Outputs: grid.pickle (incremental)
#
# To Run:
# ./create_grid.py
#
# Modification History:
#   - v1.0 : first definition, Marco Rocchetto, December 2014
#
###########################################################



#loading libraries
import numpy
import pylab
import sys
import os
import optparse
import time
import logging
import pickle
import datetime

from numpy import * #nummerical array library
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

#trying to initiate MPI parallelisation
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    MPIsize = MPI.COMM_WORLD.Get_size()
    if MPIsize > 1:
        MPIimport = True
        if MPIrank == 0:
            MPImaster = True
            MPIverbose= True
        else:
            MPImaster = False
            MPIverbose= False
    else:
        MPIverbose = True
        MPImaster  = True
        MPIimport  = False
except ImportError:
    MPIimport = False
    MPImaster = True
    MPIverbose= True

#checking for multinest library
try:
    import pymultinest
    multinest_import = True
except:
    logging.warning('Multinest library is not loaded. Multinest functionality will be disabled')
    multinest_import = False

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from tp_profile import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

# Define some constants

RSOL = 6.955e8 # solar radius in m
RJUP = 6.9911e7 # jupiter radius in m


# get arguments
parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  help='Parameter filename')
parser.add_option('-o', '--output',
                  dest="out_dir",
                  help='Output folder',
)
parser.add_option('-r', '--resolution',
                  dest="resolution",
                  default=False,
                  help='Resolution(s), separated by comma (if not specified, will be taken from parameter file)'
)
parser.add_option('-s', '--signaltonoise',
                  dest="snr",
                  default=False,
                  help='Signal to noise ratio(s), separated by comma (if not specified, will be taken from parameter file)'
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
                  action="store_true",
)
options, remainder = parser.parse_args()

if not options.param_filename or not options.out_dir:
    print parser.print_help()
    sys.exit()

# get parameters from .par file
params = parameters(options.param_filename)
out_dir = options.out_dir

# Import MPI
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.info('MPI disabled')

# get resolutions and snr (from arguments or par file)
if options.resolution:
    resolutions = [int(m) for m in options.resolution.split(',')]
elif params.grid_res:
    resolutions = params.grid_res
else:
    logging.error('Specify resolution(s) for your model spectra')
    exit()

if options.snr:
    snrs = [float(m) for m in options.snr.split(',')]
elif params.grid_snr:
    snrs = params.grid_snr
else:
    logging.error('Specify snr(s) for your model spectra')
    exit()


# check if output pickle file already exist
pickle_location = os.path.join(out_dir, 'grid.db')
if os.path.isfile(pickle_location):
    logging.error('The output pickle file already exists (%s). Delete it or rename it' % pickle_location)
    exit()

logging.info('Output grid will be saved into: %s' % pickle_location)
logging.info('Chains and plots for each snr/res set stored in: %s' % out_dir)

# initialise output pickle file and store planet/spectrum details
outdb = {
    'radius': params.planet_radius,
    'mass': params.planet_mass,
    'albedo': params.planet_albedo,
    'mu': params.planet_mu,
    'molecules_input': params.planet_molec,
    'mixing_input': params.planet_mixing,
    'temperature_input': params.planet_temp,
    'wavelength_range': (params.gen_wavemin, params.gen_wavemax),
    'resolutions': resolutions,
    'snrs': snrs,
    'model_fitting': {}
}

# loop over all resolutions and snrs
# each iteration is saved as a dictionary dictionary into outdb['model_fitting'][resoltion][snr]
for resolution in resolutions:

    if MPIrank == 0 and not resolution in outdb['model_fitting']:
        outdb['model_fitting'][resolution] = {}

    for snr in snrs:
        if MPIrank == 0:
            # set output location of plots and chains
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)
            if not os.path.isdir(os.path.join(out_dir, 'r%.1f-snr%1.f' % (resolution, snr))):
                os.mkdir(os.path.join(out_dir, 'r%.1f-snr%1.f' % (resolution, snr)))

        # update params output folder
        params.out_path = os.path.join(out_dir, 'r%.1f-snr%1.f' % (resolution, snr))
        params.out_save_plots = params.out_path


        # print some logs
        logging.info('Fitting for snr=%.1f and res=%.1f' % (snr, resolution))
        logging.info('Output chains and plots will be stored in %s' % params.out_path)
        logging.info('Generating spectrum with snr=%.1f and res=%.1f' % (snr, resolution))
        logging.info('Wavelength range: %.1f - %.1f micron' % (params.gen_wavemin, params.gen_wavemax))
        logging.info('Planet: R=%.3e R_J, M=%.3e M_J, T=%.1f, albedo=%.1f' % (params.planet_radius,
                                                                              params.planet_mass,
                                                                              params.planet_temp,
                                                                              params.planet_albedo))
        logging.info('Atmospheric composition: %s (X= %s)' % (', '.join(params.planet_molec),
                                                              ', '.join(format(x, ".2e") for x in params.planet_mixing)))

        # create spectrum
        # Planet & atmosphere properties are defined in the user provided par file
        params.gen_spec_res = resolution
        params.gen_manual_waverange = True
        dataob = data(params)
        profileob = tp_profile(dataob)

        #adding some molecules to the atmosphere @todo move to better place
        dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
        dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

        #initialising transmission radiative transfer code object
        if params.gen_type == 'transmission':
            transob = transmission(profileob)
            if params.trans_cpp:
                model = transob.cpath_integral()  # computing transmission
            else:
                model = transob.path_integral()  # computing transmission
        elif params.gen_type == 'emission':
            emisob = emission(profileob)
            model = emisob.path_integral()  # computing transmission

        # store simulated spectrum
        spectrum = np.zeros((len(dataob.specgrid),3))
        spectrum[:,0] = dataob.specgrid
        spectrum[:,1] = model
        planet_star_ratio = (params.planet_radius*RJUP)/(params.star_radius*RSOL)

        # adding errorbars to spectrum
        spec_amp = np.max(spectrum[:,1])-np.min(spectrum[:,1])
        sigma = spec_amp / float(snr)
        spectrum[:,2] = sigma

        # Fitting spectrum
        logging.info('Fitting spectrum')

        params.gen_spec_res = 1000
        params.gen_manual_waverange = False

        dataob = data(params, spectrum=spectrum)
        profileob = tp_profile(dataob)

        #adding some molecules to the atmosphere @todo move to better place
        dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
        dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

        #initialising transmission radiative transfer code instance
        if params.gen_type == 'transmission':
            transob = transmission(profileob)
        #initialising emission radiative transfer code instance
        elif params.gen_type == 'emission':
            emissob = emission(profileob)

        fitob = fitting(profileob)

        if params.gen_type == 'transmission':
            fitob.set_model(transob) #loading transmission model into fitting object
        elif params.gen_type == 'emission':
            fitob.set_model(emissob) #loading emission model into fitting object

        # use nested sampling to get posteriors
        fitob.multinest_fit() # Nested sampling fit

        # save output only from main thread
        if MPIrank == 0:

            # output
            outputob = output(fitob)
            #outputob.plot_all(save2pdf=params.out_save_plots)


           # store fitted spectrum
            fitted_spectra = []
            for spec in outputob.spec_nest:
                fitted_spectrum = np.zeros((len(spectrum[:,0]),2))
                fitted_spectrum[:,0] = spectrum[:,0]
                fitted_spectrum[:,1] = np.transpose(spec)
                fitted_spectra.append(fitted_spectrum)
                #if params.out_dump_internal: #saving models to ascii
                #    @todo needs to save output to ascii
                #    outputob.save_model(modelout=fitted_spectrum, modelsaveas=params.out_internal_name)

            if params.out_save_plots:
                outputob.plot_all(save2pdf=params.out_save_plots)

            fit_dict = {
                'observed_spectrum': spectrum, # observed spectrum (array)
                'fitted_spectra': fitted_spectra ,# fitted spectrum (array)
                'modes': fitob.NEST_modes, # mixing ratios as a function of atm layer
                'X': fitob.NEST_X_mean, # mixing ratios as a function of atm layer
                'X_std': fitob.NEST_X_std, # error on mixing ratio (1d array)
                'T': fitob.NEST_T_mean, # fitted temperature
                'T_std': fitob.NEST_T_std, # fitted T error
                'stats': fitob.NEST_stats,
                'dir_multinest': os.path.abspath(fitob.dir_multinest), # directory where nested sampling outputs are stored
                'datetime': datetime.datetime.now(),
            }

            outdb['model_fitting'][resolution][snr] = fit_dict
            logging.info('+++++++++++++ Fitting for snr=%.1f res=%.1f completed! +++++++++++++' % (resolution, snr))

        # wait for everyone to syncronise
        MPI.COMM_WORLD.Barrier()

if MPIrank == 0:
    logging.info('Fits for all snr/res comptued! Saving results to pickle file: %s' % pickle_location)
    pickle.dump(outdb, open(pickle_location, 'wb'))

logging.info('Program terminated')
