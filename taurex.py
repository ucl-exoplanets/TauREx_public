'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    TauREx execution code

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

# loading libraries
import  sys
import os
import optparse
from ConfigParser import SafeConfigParser
import logging

# setting some parameters to global
global MPIimport, MPIverbose, MPImaster, multinest_import

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

# checking for multinest library
try: 
    import pymultinest
    multinest_import = True
except:
    logging.warning('Multinest library is not loaded. Multinest functionality will be disabled')
    multinest_import = False

# loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters
from parameters import *

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

# Initialise parameters instance
params = parameters(options.param_filename)

# MPI related message
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.info('MPI disabled')

# modifying parameters object if running in cluster mode (see cluster class for docu)
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

MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

#running Tau-REx
run(params)