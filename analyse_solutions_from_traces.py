#loading libraries
import sys, os, optparse, time
import numpy as np
import pylab as pl
from ConfigParser import SafeConfigParser

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,atmosphere,data,preselector
from parameters import *
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *


parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
                  action="store_true",
)
parser.add_option('-d', '--dir',
                  dest="path",
                  default="params",
)


options, remainder = parser.parse_args()
params = parameters(options.param_filename)

if options.path is "params":
    out_path_orig = params.out_path
else:
    out_path_orig = options.path

params.out_save_plots = True

if params.gen_type == 'transmission' or params.fit_transmission:
    if os.path.isdir(os.path.join(out_path_orig, 'stage_0')):
        params.out_path = os.path.join(out_path_orig, 'stage_0')
    else:
        params.out_path = out_path_orig
            
    dataob = data(params)
    atmosphereob = atmosphere(dataob)
    forwardmodelob = transmission(atmosphereob)
    fittingob = fitting(forwardmodelob)
    if params.mcmc_run and pymc_import:
        fittingob.MCMC = True
    if params.nest_run and multinest_import:
        fittingob.NEST = True
    outputob = output(fittingob)
    if params.verbose or params.out_save_plots:
        outputob.plot_all(save2pdf=params.out_save_plots, param_labels=fittingob.fit_params_texlabels)
    outputob.save_ascii_spectra()
    
    
if params.gen_type == 'emission' or params.fit_emission:
    
    folders = ['stage_0', 'stage_1']
    for f in folders:
        dir = os.path.join(out_path_orig, f)
        if os.path.isdir(dir):
            params.out_path = dir
            dataob = data(params)
            if f is 'stage_1':
                Cov_array = np.loadtxt(os.path.join(out_path_orig, 'stage_0/tp_covariance.dat'))
                atmosphereob = atmosphere(dataob, tp_profile_type='hybrid', covariance=Cov_array)
            else:
                atmosphereob = atmosphere(dataob) 
            forwardmodelob = emission(atmosphereob)
            fittingob = fitting(forwardmodelob)
            if params.mcmc_run and pymc_import:
                fittingob.MCMC = True
            if params.nest_run and multinest_import:
                fittingob.NEST = True
            outputob = output(fittingob)
            if params.verbose or params.out_save_plots:
                outputob.plot_all(save2pdf=params.out_save_plots,
                               params_names=fittingob.fit_params_names[:fittingob.fit_X_nparams],
                               params_labels=fittingob.fit_params_texlabels[:fittingob.fit_X_nparams])
            outputob.save_ascii_spectra()
            # save and plot TP profile (plotting only if save2pdf=True)
            outputob.save_TP_profile(save2pdf=True)  #saving TP profile
            
            #deleting objectes
            dataob = None; atmosphereob = None; forwardmodelob = None; fittingob = None;
            outputob = None
            del dataob; del atmosphereob; del forwardmodelob; del fittingob; del outputob;

    
    
    