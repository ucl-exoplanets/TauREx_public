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
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    if MPIrank == 0:
        MPImaster = True
    else:
        MPImaster = False
except ImportError:
    MPIimport = False
    MPImaster = True
    

#checking for multinest library
try: 
    import pymultinest
    multinest_import = True
except:
    print 'WARNING: Multinest library is not loaded. Multinest functionality will be disabled.'
    multinest_import = False


#checking for pymc library
global pymc_import
try: 
    import pymc
    pymc_import = True
except:
    print 'WARNING: PYMC library is not loaded. MCMC functionality will be disabled.'
    pymc_import = False


#start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()

#end of profiling code

#loading classes
from classes.parameters import *
# from classes.emission import *
from classes.transmission import *
from classes.output import *
from classes.fitting import *
from classes.profile import *
from classes.data import *
from classes.preselector import *

#loading libraries
# from library.library_emission import *
from library.library_transmission import *
from library.library_general import *
from library.library_plotting import *


parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true",
)
options, remainder = parser.parse_args()

#Initialise parameters object
params = parameters(options.param_filename)
# if params.verbose:
#     print 'ARGV      :', sys.argv[1:]
#     print 'VERBOSE   :', params.verbose
#     print 'PARFILE    :', options.param_filename
#     print 'REMAINING :', remainder

#####################################################################

##check for MPI support
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.Get_size()
# print("my rank is %d"%rank)
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


#initialising data object
dataob = data(params)

#initiating preselector class
if params.pre_run:
    if params.verbose: print 'loading preprocessing'
    preob = preselector(params,dataob)
    # preob.run_preprocess(convertLinelist=False,generateSpectra=True,generatePCA=True)
    preob.run()
    if params.verbose: print preob.molselected
    params = preob.update_params()
    dataob.reset(params)


#adding some molecules to the atmosphere
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

#printing a few data object attributes
# print dataob.atmosphere['mol']
# print dataob.atmosphere['info']['mu']
# print dataob.atmosphere['mol'].keys()

#initialising TP profile object
if params.verbose: print 'loading profile'
profileob = profile(params, dataob)

#initialising transmission radiative transfer code object
if params.verbose: print 'loading transmission'
transob = transmission(params, dataob)


#initialising fitting object
if params.verbose: print 'loading fitting'
fitob = fitting(params, dataob, profileob, transob)


if params.verbose: print 'fitting data'
#fit data
if params.fit_transmission:
    fitob.downhill_fit()    #simplex downhill fit
    if params.mcmc_run and pymc_import:
#     if params.mcmc_run and pymc_import and MPImaster:
        fitob.mcmc_fit()    #MCMC fit
    if params.nest_run and multinest_import:
        fitob.multinest_fit()   #Nested sampling fit


#
outputob = output(params, dataob, fitob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
if params.verbose: outputob.plot_all()
# outputob.plot_spectrum()   #plotting data only
# outputob.plot_multinest()  #plotting multinest posteriors
# outputob.plot_mcmc()       #plotting mcmc posterios
# outputob.plot_fit()        #plotting model fits
#
outputob.save_model()       #saving models to ascii
#
#
#


#profiling code
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
if params.verbose: show()