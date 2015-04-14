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

###start of profiling code
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

AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-r', '--res',
                  dest="resolution",
                  default=1000,
)
parser.add_option('-n', '--noise',
                  dest="noise",
                  default=0,
)
parser.add_option('-e', '--error',
                  dest="error",
                  default=50,
)
parser.add_option('-T', '--T_profile',
                  dest="tp_profile",
                  default=0,
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

# set model resolution to 1000
params.gen_spec_res = 1000
params.gen_manual_waverange = True

#initialising data object
dataob = data(params)

#initialising TP profile object
atmosphereob = atmosphere(dataob)

## Apply TP profile
def movingaverage(values,window):
    weigths = np.repeat(1.0, window)/window
    smas = np.convolve(values, weigths, 'valid')
#     smas2 = np.convolve(smas[::-1],weigths,'valid')
    return smas #smas2[::-1] # as a numpy array
if int(options.tp_profile) == 1:
    logging.info('Applying custom TP profile')
    MAX_P = atmosphereob.P[0]
    MIN_P = atmosphereob.P[-1]
    smooth_window = 5 #smoothing window size as percent of total data
    Pnodes = [MAX_P,1e5, 500.0,MIN_P]
    Tnodes = [2200,2200, 1700,1700]
    TP = np.interp((np.log(atmosphereob.P[::-1])), np.log(Pnodes[::-1]), Tnodes[::-1])
    #smoothing T-P profile
    wsize = atmosphereob.nlayers*(smooth_window/100.0)
    if (wsize %2 == 0):
        wsize += 1
    TP_smooth = movingaverage(TP,wsize)
    border = np.int((len(TP) - len(TP_smooth))/2)
    atmosphereob.T = TP[::-1]
    atmosphereob.T[border:-border] = TP_smooth[::-1]

    out = np.zeros((len(atmosphereob.T),2))
    out[:,0] = atmosphereob.T
    out[:,1] = atmosphereob.P
    np.savetxt(os.path.join(params.out_path, 'TP_profile.dat'), out)
    
    figure()
    plot(atmosphereob.T, atmosphereob.P)
    pl.plot(Tnodes,Pnodes, 'x')
    pl.xlim(np.min(Tnodes)-np.min(Tnodes)*0.1,np.max(Tnodes)+np.max(Tnodes)*0.1)
    pl.yscale('log')
    pl.xlabel('Temperature')
    pl.ylabel('Pressure (Pa)')
    pl.gca().invert_yaxis()

print 'The mean molecular weight is', atmosphereob.planet_mu/AMU



#initialising transmission radiative transfer code object
if params.gen_type == 'transmission':
    forwardmodelob = transmission(atmosphereob)
elif params.gen_type == 'emission':
    forwardmodelob = emission(atmosphereob)

# bin down internal model to given resolution (default = 1000)
wavegrid, dlamb_grid = dataob.get_specgrid(R=int(options.resolution),lambda_min=params.gen_wavemin,lambda_max=params.gen_wavemax)
spec_bin_grid, spec_bin_grid_idx = dataob.get_specbingrid(wavegrid, dataob.specgrid)
model = forwardmodelob.model()
model_binned = [model[spec_bin_grid_idx == i].mean() for i in xrange(1,len(spec_bin_grid))]

out = np.zeros((len(wavegrid),3))
out[:,0] = wavegrid
out[:,1] = model_binned
out[:,2] += float(options.error) * 1e-6

if int(options.noise) == 1:
    out[:,1] += np.random.normal(0, float(options.error) * 1e-6, len(wavegrid))


#initiating output object with fitted data from fitting class
outputob = output(forwardmodel=forwardmodelob)

#plotting fits and data
outputob.plot_manual(out, save2pdf=params.out_save_plots)

if params.out_dump_internal:
    #saving models to ascii
    outputob.save_spectrum_to_file(spectrum=out, saveas=params.out_internal_name)


#end of main code
#####################################################################



###profiling code
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

# ps.print_stats()
# print s.getvalue()


#last line. displays any diagrams generated. must be run after profiling
if params.verbose:
    show()
