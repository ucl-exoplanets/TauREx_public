###########################################################
# TauREx emission execution code.     
###########################################################
#loading libraries     
import sys,os, logging
import numpy as np #nummerical array library 
import pylab as pl #science and plotting library for python

#loading classes
sys.path.append('./classes')
import emission, output, fitting, atmosphere, data, preselector
from emission import *
from output import *
from fitting import *
from atmosphere import *
from data import *
from preselector import *

import pickle

#loading libraries
sys.path.append('./library')
import library_emission as emlib 

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
#     fittingob.downhill_fit()    #simplex downhill fit
# 
    if params.mcmc_run and pymc_import:
        fittingob.mcmc_fit() # MCMC fit
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here
  
    if params.nest_run and multinest_import:
        fittingob.multinest_fit() # Nested sampling fit
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

    
    #loading output module for stage 0
    outputob = output(fittingob,out_path=os.path.join(out_path_orig, 'stage_0'))
    
    
    
    #generating TP profile covariance from previous fit
    Cov_array = emlib.generate_tp_covariance(outputob)

    #saving covariance
    if MPIimport and MPI.COMM_WORLD.Get_rank() is 0 or MPIimport is False:
        np.savetxt(os.path.join(params.out_path, 'tp_covariance.dat'), Cov_array)
    # Cov_array = np.loadtxt(os.path.join(params.out_path, 'tp_covariance.dat'))

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
        #MPI.MPI_Finalize()
        exit()
        
    #initiating output instance with fitted data from fitting class
    if params.fit_emission_stage2:
        outputob1 = output(fittingob1,out_path=os.path.join(out_path_orig, 'stage_1'))

    #plotting fits and data
    logging.info('Plotting and saving results')

    if params.verbose or params.out_save_plots:
        outputob.plot_all(save2pdf=params.out_save_plots,
                          params_names=fittingob.fit_params_names[:fittingob.fit_X_nparams])
        if params.fit_emission_stage2:
            outputob1.plot_all(save2pdf=params.out_save_plots,
                               params_names=fittingob.fit_params_names[:fittingob.fit_X_nparams])
     


    outputob.save_ascii_spectra()       #saving models to ascii
    if params.fit_emission_stage2:
        outputob1.save_ascii_spectra()       #saving models to ascii

    # save and plot TP profile (plotting only if save2pdf=True)
    outputob.save_TP_profile(save2pdf=params.out_save_plots)  #saving TP profile
    if params.fit_emission_stage2:
        outputob1.save_TP_profile(save2pdf=params.out_save_plots)


