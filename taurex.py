#! /usr/bin/python -W ignore

###########################################################
# TauRex (formerly Exonest) - Inverse spectral retrieval code using nested sampling
#
# Requirements: -python libraries: pylab, numpy, ConfigParser 
#             [these are the minimum requirements]
#
# Additional requirements: pymc
#
# Inputs: 
#
# Outputs: 
#
# To Run: 
# ./exonest.py [-p exonest.par] 
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013       
#       
###########################################################

#loading libraries     
import numpy, pylab, sys, os, optparse, time
from numpy import * #nummerical array library 
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser
import logging

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


#checking for pymc library
global pymc_import
try: 
    import pymc
    pymc_import = True
except:
    logging.warning('PYMC library is not loaded. MCMC functionality will be disabled')
    pymc_import = False


####start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()

#end of profiling code

##check for MPI support
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.Get_size()
# print "my rank is %d"%rank, ' ',MPIverbose
# 
# 
# exit()
# 
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# # size=comm.size
# print("my rank is %d"%rank)

# print 'MyRank ', MPIrank
# if MPImaster:
#     print 'MasterRank ', MPIrank
# exit()
# MPI.Init()

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


#loading parameter file parser
parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="Parfiles/exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true",
)
options, remainder = parser.parse_args()

#Initialise parameters instance
params = parameters(options.param_filename)


#####################################################################
#beginning of main code

#MPI related message
if MPIverbose and MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
elif MPIverbose:
    logging.info('MPI disabled')

# initialising data object
dataob = data(params)

#adding bulk composition to the atmosphere @todo move to better place
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

#initialising TP profile instance
profileob = tp_profile(dataob)

#initiating and running preselector instance
if params.pre_run:

    if params.fit_transmission:

        model_object = transmission(profileob)
        preob = preselector(model_object)

    elif params.fit_emission:
        model_object = emission(params,dataob,profileob) #@todo fix it as well
        preob = preselector(model_object)

    preob.run()
    logging.info('Molecule pre-selected using Marple module: %s' % preob.molselected)
    logging.info('Updating parameters object values')
    params = preob.update_params()

    #dataob.reset(updated_params) # @todo check! tp_profile needs to be updated as well

    # here we need to reinitialise the dataob and tp_profile objects:
    # @todo However, we could move everything related to the atmospheric composition to tp_profile.
    logging.info('Reinitialising data and tp_profile objects')
    dataob = data(updated_params)
    profileob = tp_profile(dataob)


#initialising transmission radiative transfer code instance
if params.fit_transmission:
    transob = transmission(profileob)

#initialising emission radiative transfer code instance
if params.fit_emission:
    emissob = emission(profileob)
#     emissob.path_integral()

# exit()

#initialising fitting object  @todo here we should actually have the transmission or emission objects as input. Then get rid of self.set_model()
fitob = fitting(profileob)

if params.fit_transmission:
    fitob.set_model(transob) #loading transmission model into fitting object
elif params.fit_emission:
    fitob.set_model(emissob) #loading emission model into fitting object

#fit data
if params.downhill_run:
    # @todo should we run the downhill fit only on the first thread? Or maybe use different starting points
    fitob.downhill_fit()    #simplex downhill fit

if params.mcmc_run and pymc_import:
    fitob.mcmc_fit() # MCMC fit
    MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

if params.nest_run and multinest_import:
    fitob.multinest_fit() # Nested sampling fit

#forcing slave processes to exit at this stage
if MPI.COMM_WORLD.Get_rank() != 0:
    exit()

#initiating output instance with fitted data from fitting class
outputob = output(fitob)

#plotting fits and data
logging.info('Plotting and saving results')

if params.verbose or params.out_save_plots:
    outputob.plot_all(save2pdf=params.out_save_plots)

# outputob.plot_spectrum()   #plotting data only
# outputob.plot_multinest()  #plotting multinest posteriors
# outputob.plot_mcmc()       #plotting mcmc posterios
# outputob.plot_fit()        #plotting model fits
#
outputob.save_model()       #saving models to ascii


#end of main code
#####################################################################



####profiling code
# pr.disable()
# 
# PROFDIR = 'Profiling/'
# if not os.path.isdir(PROFDIR):
#         os.mkdir(PROFDIR)
# 
# # s = StringIO.StringIO()
# sortby = 'cumulative'
# # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# 
# #printing to terminal
# # globalstats=pstats.Stats(pr).strip_dirs().sort_stats("cumulative")
# # globalstats.print_stats()
# 
# #redirecting output to files
# sys.stdout = open(PROFDIR+'gprofile_cum.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("cumulative")
# globalstats.print_stats()
# 
# sys.stdout = open(PROFDIR+'gprofile_ncalls.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("ncalls")
# globalstats.print_stats()
# 
# sys.stdout = open(PROFDIR+'gprofile_module.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("module")
# globalstats.print_stats()
# 
# sys.stdout = open(PROFDIR+'gprofile_tottime.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("tottime")
# globalstats.print_stats()
# 
# sys.stdout = open(PROFDIR+'gprofile_pcalls.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("pcalls")
# globalstats.print_stats()
#
# # ps.print_stats()
# # print s.getvalue()


#last line. displays any diagrams generated. must be run after profiling
if params.verbose:
    show()