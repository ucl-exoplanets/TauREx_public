'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    TauREx transmission execution code

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

import sys
import os
import numpy as np #nummerical array library 
import logging

 #loading classes
sys.path.append('./classes')

import transmission, output, fitting, atmosphere, data
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *

try:
    from mpi4py import MPI

    MPIimport = True
except ImportError:
    MPIimport = False


def run(params, options=False, gpu=False):

    # initialising data object
    dataob = data(params)
          
    #initialising TP profile instance
    atmosphereob = atmosphere(dataob)
        
    #initialising transmission radiative transfer code instance
    forwardmodelob = transmission(atmosphereob, gpu=gpu)
        
    #initialising fitting object 
    fittingob = fitting(forwardmodelob)


    #fit data
    if params.downhill_run:
        if MPIimport:
            if MPI.COMM_WORLD.Get_rank() == 0:
                fittingob.downhill_fit()  # simplex downhill fit, only on first core
        else:
            fittingob.downhill_fit()  # simplex downhill fit, only on first core



    if MPIimport:
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

    if params.mcmc_run and pymc_import:
        fittingob.mcmc_fit() # MCMC fit
        if MPIimport:
            MPI.COMM_WORLD.Barrier()


    if (not options and params.nest_run and multinest_import) \
        or (params.nest_run and multinest_import and not options.no_multinest):
        fittingob.multinest_fit() # Nested sampling fit
    elif options and (params.nest_run and multinest_import and options.no_multinest):
        fittingob.NEST = True

    # exit if the rank of MPI process is > 0 (.e. leave only master process running)
    # Useful to get just one output and not N (where N is the number of processes)
    if MPIimport:
        MPIsize = MPI.COMM_WORLD.Get_size()
        if MPI.COMM_WORLD.Get_rank() > 0:
            sys.exit()
    else:
        MPIsize = 0

    # initiating output instance with fitted data from fitting class
    # running inteprolations with nthreads = MPIsize
    outputob = output(fittingob)

    return outputob