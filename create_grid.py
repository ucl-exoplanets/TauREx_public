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

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-r', '--resolution',
                  dest="resolution",
                  default=50,
)
parser.add_option('-s', '--signaltonoise',
                  dest="snr",
                  default=10,
)
parser.add_option('-o', '--output',
                  dest="output",
                  default=False,
)
parser.add_option('-i', '--incremenetal',
                  dest="incremental_db",
                  default=False,
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
                  action="store_true",
)
options, remainder = parser.parse_args()

params = parameters(options.param_filename)

if options.output:
    params.out_path = options.output


#####################################################################
#beginning of main code

#MPI related message
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.info('MPI disabled')


resolutions = [int(m) for m in options.resolution.split(',')]
snrs = [float(m) for m in options.snr.split(',')]

for resolution in resolutions:
    for snr in snrs:

        # Fitted parameters will be stored in a dictionary stored into a pickle file.
        # This dictionary (and the pickle file) can be incrementally populated if -inc=True

        # set output pickle location.
        pickle_location = os.path.join(params.out_path, 'create_grid', 'grid.pickle')
        if options.incremental_db == 'True':
            pickle_incremental = True
        else:
            pickle_incremental = False
        logging.info('Fitted parameters will be stored in a dictionary stored into a pickle file: %s' % pickle_location)
        logging.info('Incremental saving set to %s' % pickle_location)

        # set output location of plots and chains
        if not os.path.isdir(os.path.join(params.out_path, 'create_grid')):
            os.mkdir(os.path.join(params.out_path, 'create_grid'))
        if not os.path.isdir(os.path.join(params.out_path, 'create_grid', params.planet_name)):
            os.mkdir(os.path.join(params.out_path, 'create_grid', params.planet_name))
        if not os.path.isdir(os.path.join(params.out_path, 'create_grid', params.planet_name)):
            os.mkdir(os.path.join(params.out_path, 'create_grid', params.planet_name))
        if not os.path.isdir(os.path.join(params.out_path, 'create_grid', params.planet_name, 'r%.1f-snr%1.f' % (resolution, snr))):
            os.mkdir(os.path.join(params.out_path, 'create_grid', params.planet_name, 'r%.1f-snr%1.f' % (resolution, snr)))
        params.out_path = os.path.join(params.out_path, 'create_grid', params.planet_name, 'r%.1f-snr%1.f' % (resolution, snr))

        logging.info('Output chains and plots will be stored in %s' % params.out_path)

        logging.info('Generating spectrum with SNR=%.1f and R=%i' % (snr, resolution))

        logging.info('Wavelength range: %.1f - %.1f micron' % (params.gen_wavemin, params.gen_wavemax))

        logging.info('Planet (%s): R=%.3e R_J, M=%.3e M_J, T=%.1f, albedo=%.1f' % (params.planet_name,
                                                                                   params.planet_radius,
                                                                                   params.planet_mass,
                                                                                   params.planet_temp,
                                                                                   params.planet_albedo))
        logging.info('Atmospheric composition: %s (X= %s)' % (', '.join(params.planet_molec),
                                                              ', '.join(format(x, ".2e") for x in params.planet_mixing)))

        # create spectrum
        # Planet & atmosphere properties are defined in the user provided par file

        params.gen_spec_res = resolution
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

        spectrum = np.zeros((len(dataob.specgrid),3))
        spectrum[:,0] = dataob.specgrid
        spectrum[:,1] = model

        # adding errorbars to spectrum
        planet_star_ratio = (params.planet_radius*RJUP)/(params.star_radius*RSOL)
        spec_mean_amp = np.mean(spectrum[:,1] - planet_star_ratio)
        sigma = spec_mean_amp / snr
        spectrum[:,2] = sigma

        logging.info('Fitting spectrum')

        # Fitting spectrum
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

        # wait for everyone to syncronise
        MPI.COMM_WORLD.Barrier() # @todo needed?

        # save output only from main thread
        if MPI.COMM_WORLD.Get_rank() == 0:


            # output
            outputob = output(fitob)

            outputob.plot_all(save2pdf=params.out_save_plots)

            logging.info('Saving results to pickle file: %s' % pickle_location)

            # store fitted spectrum
            fitted_spectrum = np.zeros((len(spectrum[:,0]),2))
            fitted_spectrum[:,0] = spectrum[:,0]
            fitted_spectrum[:,1] = np.transpose(outputob.spec_nest)

            # store output in dictionary, then store to pickle file
            # the pickle is a list of dictionaries

            if pickle_incremental == True:
                if os.path.isfile(pickle_location):
                    pickle_file = pickle.load(open(pickle_location, 'rb'))
                else:
                    pickle_file = []
            else:
                pickle_file = []

            planet_dict = {}
            planet_dict['name'] = params.planet_name
            planet_dict['snr'] = snr
            planet_dict['resolution'] = resolution
            planet_dict['molecules'] = params.planet_molec
            planet_dict['spectrum'] = spectrum
            planet_dict['fitted_spectrum'] = fitted_spectrum
            planet_dict['X'] = fitob.NEST_X_mean
            planet_dict['X_std'] = fitob.NEST_X_std
            planet_dict['T'] = fitob.NEST_T_mean
            planet_dict['T_std'] = fitob.NEST_T_std
            planet_dict['NEST_OUT'] = fitob.NEST_FITDATA
            planet_dict['datetime'] = datetime.datetime.now()

            pickle_file.append(planet_dict)

            pickle.dump(pickle_file, open(pickle_location, 'wb'))