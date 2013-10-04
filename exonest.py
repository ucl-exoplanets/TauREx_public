#! /usr/bin/python -W ignore

###########################################################
# ExoNest - Inverse spectral retrieval code using nested sampling
#
# Requirements: -python libraries: pylab, numpy, ConfigParser 
#             [these are the minimum requirements]
#
# Additional requirements: pymc
#
# Inputs: 
#
# Outputs: 
#
# To Run: 
# ./obs_pipeline.py [-p exonest.par] [-v] 
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013       
#       
###########################################################

#loading libraries     
import numpy, pylab, sys, os, optparse, time
from numpy import * #nummerical array library 
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

starttime = time.clock()

#loading classes
from classes.parameters import *
from classes.emission import *
from classes.transmission import *
from classes.fitting import *
from classes.profile import *
from classes.data import *

#loading libraries
from library.library_emission import *
from library.library_transmission import *
from library.library_general import *


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
if params.verbose == True:
    print 'ARGV      :', sys.argv[1:]
    print 'VERBOSE   :', params.verbose
    print 'PARFILE    :', options.param_filename
    print 'REMAINING :', remainder

#####################################################################


#initialising data object
dataob = data(params)

#adding some molecules to the atmosphere
dataob.add_molecule('H2',2.0,2.0e-9,1.0001384,0.85)
dataob.add_molecule('He',4.0,1.0e-9,1.0000350,0.15)

print dataob.atmosphere['mol']

print dataob.atmosphere['info']['mu']
print dataob.atmosphere['mol'].keys()

profileob = profile(params,dataob)



transob = transmission(params,dataob)

# if params.trans_cpp == True:
#     absorption = transob.cpath_integral(rho=profileob.get_rho())
# else:
#     absorption = transob.path_integral()


fitob = fitting(params,dataob,profileob,transob)

# ble = fitob.downhill_fit()
ble = fitob.mcmc_fit()

print ble.stats()

absorption = transob.cpath_integral(rho=profileob.get_rho(T=ble[0]),X=ble[1:].reshape(fitob.Pshape[0],fitob.Pshape[1]-1))


# print absorption

figure(1)
errorbar(dataob.spectrum[:,0],dataob.spectrum[:,1],dataob.spectrum[:,2])
plot(dataob.spectrum[:,0],transpose(absorption),c='r')


show()