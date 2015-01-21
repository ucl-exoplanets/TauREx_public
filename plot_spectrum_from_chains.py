#! /usr/bin/python -W ignore

###########################################################
# plot_spectrum_from_chains: plots spectra from MCMC and/or
# Nested Sampling runs. It is based on the create_spectrum.py
# and plot_chains.py codes and may be integrated in the future. 
#
# Inputs: minimal or normal parameter file, chains directory
#
# Outputs: spectrum.dat
#
# To Run: 
# ./create_spectrum.py [-p exonest.par] [-n plot nested only] 
#       [=m plot mcmc only] [-a plot everything]
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, June 2014      
#       
###########################################################

#loading libraries     
import numpy, pylab, sys, os, optparse, time
from numpy import * #nummerical array library 
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

####start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()


#checking for multinest library
global multinest_import 
try: 
    with os.devnull as sys.stdout: 
        import pymultinest
        multinest_import = True
except:
    multinest_import = False

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,atmosphere,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from atmosphere import *
from data import *
from preselector import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *



parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="Parfiles/taurex.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true",
)
parser.add_option('-n', '--nested',
                  dest="plot_nest",
                  default=False,
                  action="store_true",
)
parser.add_option('-m', '--mcmc',
                  dest="plot_mcmc",
                  default=False,
                  action="store_true",
)
parser.add_option('-a', '--all',
                  dest="plot_all",
                  default=False,
                  action="store_true",
)
parser.add_option('-i', '--inputdir',
                  dest="chaindir",
                  default='Output/',
                  action="store",
)
parser.add_option('-o', '--outdir',
                  dest="outdir",
                  default='./Output',
                  action="store",
)

options, remainder = parser.parse_args()

print options.param_filename

#Initialise parameters object
params = parameters(options.param_filename)


#####################################################################
#getting mixing ratios and temperatures from chain directory
#tested for transmission. NOT WORKING FOR EMISSION AS PARAMETERS ARE HARDCODED!!!
#FIX IT MARCO FIX IT!

chaindir = options.chaindir
plot_mcmc = options.plot_mcmc 
plot_nest = options.plot_nest
plot_all = False

params.out_path          = options.outdir
params.out_dump_internal = True
params.pre_run           = False
params.fit_transmission  = False
params.fit_emission      = False
params.mcmc_run          = False
params.nest_run          = False


fit_params = np.loadtxt(chaindir+'parameters.dat',dtype='string')
# print fit_params

params.planet_mixing = []

#setting list of molecules
params.planet_molec = []
for mol in fit_params[:-1]:
    params.planet_molec.append(mol)


#loading MCMC state file 
mcmcdir = chaindir+'MCMC/'
if plot_mcmc: 
    params.out_internal_name = 'mcmc_' +params.out_internal_name
    with open(mcmcdir+'thread_0/state.txt','r') as statefile:
                mcmcstate = eval(statefile.read())

    params_mcmc = mcmcstate['stochastics'].keys()
    params_mcmc[1:] = np.sort(params_mcmc[1:])
    
    print params_mcmc
    
    #setting planetary temperature
    params.planet_temp = mcmcstate['stochastics'][params_mcmc[-1]]
#     Ampshift = mcmcstate['stochastics'][params_mcmc[1]]
    
    #setting mixing ratios
    for mix in params_mcmc[2:]:
        params.planet_mixing.append(float(mcmcstate['stochastics'][mix]))

    
nestdir = chaindir+'multinest/'
if plot_nest:
    params.out_internal_name = 'nested_' +params.out_internal_name
    import pymultinest as pymn
    #reading in Nested data
    nest_raw = pymn.Analyzer(n_params=len(fit_params),outputfiles_basename=nestdir+'1-')
    neststate = nest_raw.get_stats()
    
    #setting planetary temperature
    params.planet_temp = neststate['modes'][0]['mean'][-1]
    
#     Ampshift = neststate['modes'][0]['mean'][1]
    
    print neststate['modes'][0]['mean'][:-1]
    
    #setting mixing ratios
    for mix in neststate['modes'][0]['mean'][:-1]:
        params.planet_mixing.append(mix)
    
print params.planet_temp
print params.planet_molec
print params.planet_mixing


#converting to array
params.planet_mixing =  np.array(params.planet_mixing)  
params.planet_molec  =  np.array(params.planet_molec) 



#####################################################################
#initialising data object
dataob = data(params)

#initialising TP profile object
if params.verbose: print 'loading atmosphere'
atmob = atmosphere(dataob)

#adding some molecules to the atmosphere
# atmob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
# atmob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)



#initialising transmission radiative transfer code object
if params.gen_type == 'transmission':
    if params.verbose: print 'loading transmission'
    transob = transmission(atmob)

    if params.trans_cpp:
        MODEL = transob.cpath_integral()  # computing transmission
    else:
        MODEL = transob.path_integral()  # computing transmission
        
#initialising transmission radiative transfer code object
if params.gen_type == 'emission':
    if params.verbose: print 'loading emission'
    emisob = emission(atmob)
    
    MODEL = emisob.path_integral()  # computing transmission        
    
# # 
OUT = np.zeros((len(dataob.specgrid),2))
OUT[:,0] = dataob.specgrid
OUT[:,1] = MODEL
# OUT[:,2] += 1e-5 #adding errorbars. can be commented

if params.gen_type == 'emission':
    outputob = output(params=params, data=dataob,forwardmodel=emisob) #initiating output object with fitted data from fitting class
if params.gen_type == 'transmission':
    outputob = output(params=params, data=dataob,forwardmodel=transob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
outputob.plot_manual(OUT,save2pdf=params.out_save_plots)   #plotting data only
# outputob.plot_all(save2pdf=params.out_save_plots)
# outputob.plot_fit()        #plotting model fits

if params.out_dump_internal:
    outputob.save_model(modelout=OUT, modelsaveas=params.out_internal_name)       #saving models to ascii



#
if params.verbose: show()