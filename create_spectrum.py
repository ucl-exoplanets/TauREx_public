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
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=True,
                  action="store_true",
)
options, remainder = parser.parse_args()

#####################################################################

#initialising parameters object
params = parameters(options.param_filename)

#initialising data object
dataob = data(params)

#initialising TP profile object
atmosphereob = atmosphere(dataob)

#initialising transmission radiative transfer code object
if params.gen_type == 'transmission':
    forwardmodelob = transmission(atmosphereob)
elif params.gen_type == 'emission':
    forwardmodelob = emission(atmosphereob)

out = np.zeros((len(dataob.specgrid),2))
out[:,0] = dataob.specgrid
out[:,1] = forwardmodelob.model()
# out[:,2] += 1e-5 #adding errorbars. can be commented

#initiating output object with fitted data from fitting class
outputob = output(forwardmodel=forwardmodelob)


#plotting fits and data
outputob.plot_manual(out, save2pdf=params.out_save_plots)

if params.out_dump_internal:
    #saving models to ascii
    outputob.save_model(modelout=out, modelsaveas=params.out_internal_name)


#end of main code
#####################################################################



####profiling code
# pr.disable()
#  
# PROFDIR = 'Profiling/'
# if not os.path.isdir(PROFDIR):
#         os.mkdir(PROFDIR)
#  
# # s = StringIO.StringIO()
# sortby = 'cumulative'
# # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#  
# #printing to terminal
# # globalstats=pstats.Stats(pr).strip_dirs().sort_stats("cumulative")
# # globalstats.print_stats()
#  
# #redirecting output to files
# sys.stdout = open(PROFDIR+'gprofile_cum.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("cumulative")
# globalstats.print_stats()
#  
# sys.stdout = open(PROFDIR+'gprofile_ncalls.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("ncalls")
# globalstats.print_stats()
#  
# sys.stdout = open(PROFDIR+'gprofile_module.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("module")
# globalstats.print_stats()
#  
# sys.stdout = open(PROFDIR+'gprofile_tottime.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("tottime")
# globalstats.print_stats()
#  
# sys.stdout = open(PROFDIR+'gprofile_pcalls.profile','wb')
# globalstats=pstats.Stats(pr).strip_dirs().sort_stats("pcalls")
# globalstats.print_stats()
# 
# # ps.print_stats()
# # print s.getvalue()


#last line. displays any diagrams generated. must be run after profiling
if params.verbose: show()