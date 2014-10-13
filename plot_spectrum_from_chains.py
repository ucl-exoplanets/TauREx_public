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

#checking for multinest library
global multinest_import 
try: 
    with os.devnull as sys.stdout: 
        import pymultinest
        multinest_import = True
except:
    multinest_import = False

#loading classes
from classes.parameters import *
# from classes.emission import *
from classes.transmission import *
from classes.output import *
from classes.profile import *
from classes.data import *

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
                  default='chains/',
                  action="store",
)
parser.add_option('-o', '--outdir',
                  dest="outdir",
                  default='./',
                  action="store",
)

options, remainder = parser.parse_args()

print options.param_filename

#Initialise parameters object
params = parameters(options.param_filename)


#####################################################################


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
for mol in fit_params[1:]:
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
    params.planet_temp = mcmcstate['stochastics'][params_mcmc[0]]
#     Ampshift = mcmcstate['stochastics'][params_mcmc[1]]
    
    #setting mixing ratios
    for mix in params_mcmc[2:]:
        params.planet_mixing.append(float(mcmcstate['stochastics'][mix]))

    
    
if plot_nest:
    params.out_internal_name = 'nested_' +params.out_internal_name
    import pymultinest as pymn
    #reading in Nested data
    nest_raw = pymn.Analyzer(n_params=len(fit_params),outputfiles_basename=chaindir+'1-')
    neststate = nest_raw.get_stats()
    
    #setting planetary temperature
    params.planet_temp = neststate['modes'][0]['mean'][0]
    
#     Ampshift = neststate['modes'][0]['mean'][1]
    
    print neststate['modes'][0]['mean'][1:]
    
    #setting mixing ratios
    for mix in neststate['modes'][0]['mean'][1:]:
        params.planet_mixing.append(mix)
    
print params.planet_temp
print params.planet_molec
print params.planet_mixing

#converting to array
params.planet_mixing =  np.array(params.planet_mixing)  
params.planet_molec  =  np.array(params.planet_molec) 


#initialising data object
dataob = data(params)

#adding some molecules to the atmosphere
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)


#initialising TP profile object
if params.verbose: print 'loading profile'
profileob = profile(params, dataob)

#initialising transmission radiative transfer code object
if params.verbose: print 'loading transmission'
transob = transmission(params, dataob)


########
###example of how to manually reading in ABS file and computing transmission spectrum
# dataob.set_ABSfile(path='/Users/ingowaldmann/UCLlocal/REPOS/exonestpy/Input/crosssections/',
#                    filelist=['1H2-16O_300-29995_600K_1.000000.sigma.abs',
#                              '12C-1H4_300-11999_600K_1.000000.sigma.abs'],interpolate=True)
# transob.reset(dataob) #resets transob to reflect changes in dataob
########



# # #manually setting mixing ratio and T-P profile
# X_in   = zeros((5,profileob.nlayers))
# # 
#gj3470 composition
# X_in[0,:]  += 5.0e-3
# X_in[1,:]  += 4.0e-3
# X_in[2,:]  += 2e-3
# X_in[3,:]  += 7e-4 #2e-7
# X_in[4,:]  += 2e-4

#carbon poor 1solar c/o=0.5
# X_in[0,:]  += 4.0e-4
# X_in[1,:]  += 2.0e-8
# X_in[2,:]  += 4.0e-4
# X_in[3,:]  += 2.0e-7 #2e-7
# X_in[4,:]  += 1.0e-5
# 
# #carbon poor 10x solar c/o=0.5
# X_in[0,:]  += 4.0e-4
# X_in[1,:]  += 1.0e-7
# X_in[2,:]  += 4.0e-4
# X_in[3,:]  += 7.0e-5 #2e-7
# X_in[4,:]  += 1.0e-5


# rho_in = profileob.get_rho(T=params.planet_temp)




if params.trans_cpp:
    MODEL = transob.cpath_integral()  # computing transmission
else:
    MODEL = transob.path_integral()  # computing transmission


#shifitng MODEL by Ampshift
# MODEL += Ampshift
    
# # 
OUT = np.zeros((len(dataob.specgrid),3))
OUT[:,0] = dataob.specgrid
OUT[:,1] = MODEL

# OUT[:,2] += 5e-5 #adding errorbars. can be commented


outputob = output(params, dataob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
# outputob.plot_manual(OUT)   #plotting model only

# if params.out_dump_internal:
outputob.save_model(modelout=OUT, modelsaveas=params.out_internal_name)       #saving models to ascii
#
if params.verbose: show()