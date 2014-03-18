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
# from classes.emission import *
from classes.transmission import *
from classes.output import *
from classes.fitting import *
from classes.profile import *
from classes.data import *
from classes.preselector import *

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
if params.verbose == True:
    print 'ARGV      :', sys.argv[1:]
    print 'VERBOSE   :', params.verbose
    print 'PARFILE    :', options.param_filename
    print 'REMAINING :', remainder

#####################################################################

#initiating preselector class
# preob = preselector(params)


# exit()

#initialising data object
dataob = data(params)

#adding some molecules to the atmosphere
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

#printing a few data object attributes
print dataob.atmosphere['mol']
print dataob.atmosphere['info']['mu']
print dataob.atmosphere['mol'].keys()

#initialising TP profile object
profileob = profile(params, dataob)

#initialising transmission radiative transfer code object
transob = transmission(params, dataob)

# #manually reading in ABS file and computing transmission spectrum
dataob.set_ABSfile(path='/Users/ingowaldmann/UCLlocal/REPOS/exonest/exonestpy/test-code/crosssections/',
                   filelist=['1H2-16O_0-29999_600K_0.010000.sigma.abs'])
transob.reset(dataob)
#
print 'ble'
MODEL = transob.cpath_integral()
print MODEL

print size(MODEL)
print size(dataob.wavegrid)
print size(dataob.spectrum)


figure()
plot(MODEL)
show()

# #initialising fitting object
# fitob = fitting(params, dataob, profileob, transob)
#
# #fit data
# # fitob.downhill_fit()    #simplex downhill fit
# # fitob.mcmc_fit()        #MCMC fit
# fitob.multinest_fit()   #Nested sampling fit
#
# #manually call transmission spectrum code
# # absorption = transob.cpath_integral(rho=profileob.get_rho(T=fitob.MCMC_T_mean),X=fitob.MCMC_X_mean)
#
#
#
# outputob = output(params, dataob, fitob) #initiating output object with fitted data from fitting class
#
# #plotting fits and data
# outputob.plot_all()
# # outputob.plot_spectrum()   #plotting data only
# # outputob.plot_multinest()  #plotting multinest posteriors
# # outputob.plot_mcmc()       #plotting mcmc posterios
# # outputob.plot_fit()        #plotting model fits
#
# # outputob.save_model()       #saving models to ascii
#
#
#
# show()