###########################################################
# TauREx transmission execution code.     
###########################################################
#loading libraries     
import pylab,sys
import numpy as np #nummerical array library 
import pylab as pl #science and plotting library for python
import logging

#loading classes
sys.path.append('./classes')
import transmission, output, fitting, atmosphere, data, preselector
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *
from preselector import *

   
def run(params):    
    # initialising data object
    dataob = data(params)
          
    #initialising TP profile instance
    atmosphereob = atmosphere(dataob)
        
    #initiating and running preselector instance
    if params.pre_run:
        
        model_object = transmission(atmosphereob)
        preob = preselector(model_object)
        
        preob.run()
        logging.info('Molecule pre-selected using Marple module: %s' % preob.molselected)
        logging.info('Updating parameters object values')
        updated_params = preob.update_params()
        
        logging.info('Reinitialising data and atmosphereob objects')
        dataob = data(updated_params)
            
        atmosphereob = atmosphere(dataob)
        
        
    #initialising transmission radiative transfer code instance
    forwardmodelob = transmission(atmosphereob)
        
    #initialising fitting object 
    fittingob = fitting(forwardmodelob)
        
    #fit data
    if params.downhill_run:
        # @todo should we run the downhill fit only on the first thread? Or maybe use different starting points
        fittingob.downhill_fit()    #simplex downhill fit
        
    if params.mcmc_run and pymc_import:
        fittingob.mcmc_fit() # MCMC fit
        MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here
        
    if params.nest_run and multinest_import:
        fittingob.multinest_fit() # Nested sampling fit
        
    #forcing slave processes to exit at this stage
    if MPIimport and MPI.COMM_WORLD.Get_rank() != 0:
        #MPI.MPI_Finalize()
        exit()
        
    #initiating output instance with fitted data from fitting class
    outputob = output(fittingob)
        
    #plotting fits and data
    logging.info('Plotting and saving results')
        
    if params.verbose or params.out_save_plots:
        outputob.plot_all(save2pdf=params.out_save_plots)
        
    # outputob.plot_spectrum()   #plotting data only
    # outputob.plot_multinest()  #plotting multinest posteriors
    # outputob.plot_mcmc()       #plotting mcmc posterios
    # outputob.plot_fit()        #plotting model fits
    #
    outputob.save_ascii_spectra()       #saving models to ascii