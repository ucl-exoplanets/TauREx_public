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
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    if MPIrank == 0:
        MPImaster = True
    else:
        MPImaster = False
except ImportError:
    MPIimport = False




#start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()

#end of profiling code

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
# if params.verbose:
#     print 'ARGV      :', sys.argv[1:]
#     print 'VERBOSE   :', params.verbose
#     print 'PARFILE    :', options.param_filename
#     print 'REMAINING :', remainder

#####################################################################

#
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.Get_size()
# print("my rank is %d"%rank)
#
# exit()

# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.size
# print("my rank is %d"%rank)

# MPI.Init()

#initialising data object


dataob = data(params)


#initiating preselector class
if params.pre_run:
    print 'loading preprocessing'
    preob = preselector(params,dataob)
    # preob.run_preprocess(convertLinelist=False,generateSpectra=True,generatePCA=True)
    preob.run()
    print preob.molselected
    params = preob.update_params()
    dataob.reset(params)


#adding some molecules to the atmosphere
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

#printing a few data object attributes
# print dataob.atmosphere['mol']
# print dataob.atmosphere['info']['mu']
# print dataob.atmosphere['mol'].keys()

#initialising TP profile object
print 'loading profile'
profileob = profile(params, dataob)

#initialising transmission radiative transfer code object
print 'loading transmission'
transob = transmission(params, dataob)


# figure(200)
# plot(transpose(dataob.sigma_array))
# show()
# exit()

########
# ###example of how to manually reading in ABS file and computing transmission spectrum
# 
# 
# dataob.set_ABSfile(path='/Users/ingowaldmann/UCLlocal/REPOS/exonestpy/test-code/crosssections/',
#                    filelist=['1H2-16O_300-29995_1400K_1.000000.sigma.abs',
#                              '12C-1H4_300-11999_1400K_1.000000.sigma.abs',
#                              '12C-16O_300-8499_1400K_1.000000.sigma.abs',
#                              '14N-1H3_300-19999_1400K_1.000000.sigma.abs',
#                              '12C-16O2_300-12999_1400K_1.000000.sigma.abs'],interpolate=True,temperature=1400)
# # # 
# # # 
# # # dataob.set_ABSfile(path='/Users/ingowaldmann/UCLlocal/REPOS/exonestpy/test-code/crosssections/',
# # #                    filelist=['12C-16O2_300-12999_600K_1.000000.sigma.abs'],interpolate=True)
# # # 
# # 
# transob.reset(dataob) #resets transob to reflect changes in dataob
# # 
# # # # #
# # # # figure(200)
# # # # plot(dataob.specgrid,dataob.sigma_array[0,:])
# # # # plot(dataob.specgrid,dataob.sigma_array[1,:])
# # # # plot(dataob.specgrid,dataob.sigma_array[2,:])
# # # # # plot(dataob.specgrid,dataob.sigma_array[3,:])
# # # # show()
# # # # exit()
# # # # #
# # # # #
# # # #manually setting mixing ratio and T-P profile
# X_in   = zeros((5,profileob.nlayers))
# # # print np.shape(X_in)
# # # X_in[0,:]  += 1.0e-03
# # # # X_in[1,:]  += 2.54462702e-08
# # # # X_in[2,:]  += 1.14925746e-03
# # # # X_in[3,:]  += 9.96460970e-05
# # # 
# # # # X_in[0,:]  += 0.0011647246764488776
# # # # X_in[1,:]  += 2.892344980711778e-08
# # # # X_in[2,:]  += 0.0011969760859621417
# # # # X_in[3,:]  += 0.00010305931257479826
# # # 
# #gj3470 composition
# # X_in[0,:]  += 5.0e-3
# # X_in[1,:]  += 4.0e-3
# # X_in[2,:]  += 2e-3
# # X_in[3,:]  += 7e-4 #2e-7
# # X_in[4,:]  += 2e-4
# 
# #carbon poor 1solar c/o=0.5
# # X_in[0,:]  += 4.0e-4
# # X_in[1,:]  += 2.0e-8
# # X_in[2,:]  += 4.0e-4
# # X_in[3,:]  += 2.0e-7 #2e-7
# # X_in[4,:]  += 1.0e-5
# 
# #carbon poor 10x solar c/o=0.5
# X_in[0,:]  += 4.0e-4
# X_in[1,:]  += 1.0e-7
# X_in[2,:]  += 4.0e-4
# X_in[3,:]  += 7.0e-5 #2e-7
# X_in[4,:]  += 1.0e-5
# 
# # # 
# # # # 1.11622920e-03   2.54462702e-08   1.14925746e-03
# # # #    9.96460970e-05
# # # #
# # # # [1143.2703550954802, 0.0011647246764488776, 2.892344980711778e-08, 0.0011969760859621417, 0.00010305931257479826]
# #   
# rho_in = profileob.get_rho(T=1400)
# MODEL = transob.cpath_integral(rho=rho_in,X=X_in,temperature=1400)  # computing transmission
# # # 
# OUT = np.zeros((len(dataob.specgrid),3))
# OUT[:,0] = dataob.specgrid
# OUT[:,1] = MODEL
# OUT[:,2] += 5e-5
# # np.savetxt('testspec2.txt',OUT)
# # 
# figure()
# # plot(dataob.spectrum[:,0],dataob.spectrum[:,1],'g')
# errorbar(OUT[:,0],OUT[:,1],OUT[:,2],color=[0.7,0.7,0.7])
# plot(OUT[:,0],OUT[:,1],'b')
# xscale('log')
# show()
# # 
# # 
# exit()

#########


#initialising fitting object
print 'loading fitting'
fitob = fitting(params, dataob, profileob, transob)
#
print 'fitting data'
#fit data


fitob.downhill_fit()    #simplex downhill fit
# fitob.mcmc_fit()        #MCMC fit
# fitob.multinest_fit()   #Nested sampling fit

#
#manually call transmission spectrum code
# absorption = transob.cpath_integral(rho=profileob.get_rho(T=fitob.MCMC_T_mean),X=fitob.MCMC_X_mean)
#
#
#
outputob = output(params, dataob, fitob) #initiating output object with fitted data from fitting class
#
#plotting fits and data
outputob.plot_all()
# outputob.plot_spectrum()   #plotting data only
# outputob.plot_multinest()  #plotting multinest posteriors
# outputob.plot_mcmc()       #plotting mcmc posterios
# outputob.plot_fit()        #plotting model fits
#
outputob.save_model()       #saving models to ascii
#
#
#


#profiling code
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
show()