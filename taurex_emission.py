'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    TauREx emission execution code

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''
# loading libraries
import os
import sys
import logging
import numpy as np

#loading classes
sys.path.append('./classes')
import emission, output, fitting, atmosphere, data
from emission import *
from output import *
from fitting import *
from atmosphere import *
from data import *

#loading libraries
sys.path.append('./library')

from library_emission import *

def run(params):

    out_path_orig = params.out_path

    ###############################
    # STAGE 1
    ###############################

    # set output directory of stage 1
    params.out_path = os.path.join(out_path_orig, 'stage_0')

    # initialising data object
    dataob = data(params)

    #initialising TP profile instance
    atmosphereob = atmosphere(dataob)

    #initialising emission radiative transfer code instance
    forwardmodelob = emission(atmosphereob)

    #initialising fitting object
    fittingob = fitting(forwardmodelob)
    
    #fit data for stage 1
    if params.downhill_run:
        fittingob.downhill_fit()    #simplex downhill fit

    if params.mcmc_run and pymc_import:
        fittingob.mcmc_fit() # MCMC fit
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here
  
    if params.nest_run and multinest_import:
        fittingob.multinest_fit() # Nested sampling fit
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

    if MPI.COMM_WORLD.Get_rank() == 0:
        outputob = output(fittingob, out_path=os.path.join(out_path_orig, 'stage_0'))

    exit()

    # todo fix stage 2

    # generating TP profile covariance from previous fit
    Cov_array = generate_tp_covariance(outputob)

    # saving covariance
    if MPIimport and MPI.COMM_WORLD.Get_rank() is 0 or MPIimport is False:
        np.savetxt(os.path.join(params.out_path, 'tp_covariance.dat'), Cov_array)

    ###############################
    # STAGE 2
    ###############################
    
    if params.fit_emission_stage2:
        # set output directory of stage 2
        params.out_path = os.path.join(out_path_orig, 'stage_1')
    
        #setting up objects for stage 2 fitting
        dataob1 = data(params) 
    
        #setting stage 2 atmosphere object
        atmosphereob1 = atmosphere(dataob1, tp_profile_type='hybrid', covariance=Cov_array)
                
        #setting stage 2 forward model 
        forwardmodelob1 = emission(atmosphereob1)
            
        #setting stage 2 fitting object
        fittingob1 = fitting(forwardmodelob1)
    
        # #running stage 2 fit
        if params.downhill_run:
            fittingob1.downhill_fit()    #simplex downhill fit
        
        if params.mcmc_run and pymc_import:
            fittingob1.mcmc_fit() # MCMC fit
            MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here
    
        if params.nest_run and multinest_import:
            fittingob1.multinest_fit() # Nested sampling fit
            MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

    ###############
    #finished fitting. Post fitting analysis from here 
    
    #forcing slave processes to exit at this stage
    if MPIimport and MPI.COMM_WORLD.Get_rank() != 0:
        exit()
        
    #initiating output instance with fitted data from fitting class
    if params.fit_emission_stage2:
        outputob1 = output(fittingob1,out_path=os.path.join(out_path_orig, 'stage_1'))