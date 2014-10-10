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
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *
from emission import *
from transmission import *
from output import *
from fitting import *
from tp_profile import *
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
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
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
profileob = tp_profile(params, dataob)

#initialising transmission radiative transfer code object
if params.gen_type == 'transmission':
    if params.verbose: print 'loading transmission'
    transob = transmission(params, dataob,profileob)

    if params.trans_cpp:
        MODEL = transob.cpath_integral()  # computing transmission
    else:
        MODEL = transob.path_integral()  # computing transmission
        
#initialising transmission radiative transfer code object
if params.gen_type == 'emission':
    if params.verbose: print 'loading emission'
    emisob = emission(params, dataob,profileob)
    
    MODEL = emisob.path_integral()  # computing transmission        
    
# # 
OUT = np.zeros((len(dataob.specgrid),2))
OUT[:,0] = dataob.specgrid
OUT[:,1] = MODEL
# OUT[:,2] += 5e-5 #adding errorbars. can be commented


outputob = output(params, dataob,emisob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
outputob.plot_manual(OUT,save2pdf=params.out_save_plots)   #plotting data only

if params.out_dump_internal:
    outputob.save_model(modelout=OUT, modelsaveas=params.out_internal_name)       #saving models to ascii
#
if params.verbose: show()