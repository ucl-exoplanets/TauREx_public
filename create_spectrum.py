#! /usr/bin/python -W ignore

###########################################################
# create_spectrum - running minimal version of TauRex to 
# create transmission (and later emission) spectra from 
# parameter file values. 
#
# Requirements: -python libraries: pylab, numpy, ConfigParser 
#             [these are the minimum requirements]
#
# Additional requirements: pymc
#
# Inputs: minimal or normal parameter file
#
# Outputs: spectrum.dat
#
# To Run: 
# ./create_spectrum.py [-p exonest.par] 
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
options, remainder = parser.parse_args()

#Initialise parameters object
params = parameters(options.param_filename)


#####################################################################

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
transob = transmission(params, dataob,profileob)


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
    
# # 
OUT = np.zeros((len(dataob.specgrid),3))
OUT[:,0] = dataob.specgrid
OUT[:,1] = MODEL
# OUT[:,2] += 5e-5 #adding errorbars. can be commented


outputob = output(params, dataob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
outputob.plot_manual(OUT)   #plotting data only

if params.out_dump_internal:
    outputob.save_model(modelout=OUT, modelsaveas=params.out_internal_name)       #saving models to ascii
#
if params.verbose: show()