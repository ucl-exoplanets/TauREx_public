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
# ./taurex.py [-p taurex.par] 
#      
#       
###########################################################

#loading libraries     
import  sys, os, optparse, time,pylab
# from numpy import * #nummerical array library 
# from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser
import logging


####start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()

#setting some parameters to global
global MPIimport,MPIverbose,MPImaster,multinest_import,pymc_import

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
try: 
    import pymc
    pymc_import = True
except:
    logging.warning('PYMC library is not loaded. MCMC functionality will be disabled')
    pymc_import = False


# #loading classes
sys.path.append('./classes')
sys.path.append('./library')
# 
import parameters
from parameters import *
from library_general import house_keeping

#loading parameter file parser
parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="Parfiles/exonest.par"
                  )
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true"
                  )
                  
parser.add_option('-c', '--cluster',
                  dest="cluster_dictionary",
                  default="None",
                  type = "string",
                  action="store"         
                  )
parser.add_option('-i', '--cluster_index',
                  dest="cluster_procid",
                  default="None",
                  type = "string",
                  action="store"         
                  )

options, remainder = parser.parse_args()

#Initialise parameters instance
params = parameters(options.param_filename)

#MPI related message
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.info('MPI disabled')

#modifying parameters object if running in cluster mode (see cluster class for docu)
if options.cluster_dictionary is not "None":
    from cluster import cluster
    c = cluster()
    c.read_dict(dict_name=options.cluster_dictionary)
    params = c.modify_params(params,options.cluster_procid)

if params.gen_type == 'transmission' or params.fit_transmission:
    from taurex_transmission import run
elif params.gen_type == 'emission' or params.fit_emission:
    from taurex_emission import run
else:
    logging.error('Forward model selected is ambiguous')
    logging.info('Check \'type\' and \'fit_emission\', \'fit_transmission\' parameters')
    logging.info('PS: you suck at this... ')
    exit()


#running Tau-REx
run(params)


#####################################################################
#launches external housekeeping script. E.g. useful to transfer data 
#from Scratch to home 
    
if params.clean_run:
    house_keeping(params,options)

#####################################################################
#### end of profiling code
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
    pylab.show()
