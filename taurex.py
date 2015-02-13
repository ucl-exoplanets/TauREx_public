#! /usr/bin/python -W ignore

###########################################################
# TauRex (formerly Exonest) - Inverse spectral retrieval code using nested sampling
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
# ./exonest.py [-p exonest.par] 
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
import logging

#trying to initiate MPI parallelisation
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    MPIsize = MPI.COMM_WORLD.Get_size()
    if MPIsize > 1:
        MPIimport = True
        if MPIrank == 0:
            MPImaster = True
            MPIverbose= True
        else:
            MPImaster = False
            MPIverbose= False
    else:
        MPIverbose = True
        MPImaster  = True
        MPIimport  = False
except ImportError:
    MPIimport = False
    MPImaster = True
    MPIverbose= True
    


#checking for multinest library
try: 
    import pymultinest
    multinest_import = True
except:
    logging.warning('Multinest library is not loaded. Multinest functionality will be disabled')
    multinest_import = False


#checking for pymc library
global pymc_import
try: 
    import pymc
    pymc_import = True
except:
    logging.warning('PYMC library is not loaded. MCMC functionality will be disabled')
    pymc_import = False


####start of profiling code
# import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()
# starttime = time.clock()

#end of profiling code

##check for MPI support
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# size=comm.Get_size()
# print "my rank is %d"%rank, ' ',MPIverbose
# 
# 
# exit()
# 
# comm = MPI.COMM_WORLD
# rank=comm.Get_rank()
# # size=comm.size
# print("my rank is %d"%rank)

# print 'MyRank ', MPIrank
# if MPImaster:
#     print 'MasterRank ', MPIrank
# exit()
# MPI.Init()

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters, emission, transmission, output, fitting, atmosphere, data, preselector



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


#loading parameter file parser
parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="Parfiles/exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true",
)
options, remainder = parser.parse_args()

#Initialise parameters instance
params = parameters(options.param_filename)


#####################################################################
#beginning of main code

#MPI related message
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.info('MPI disabled')

# initialising data object
dataob = data(params)

#adding bulk composition to the atmosphere @todo move to better place
#dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
#dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)

#initialising TP profile instance
atmosphereob = atmosphere(dataob)

#initiating and running preselector instance
if params.pre_run:

    if params.fit_transmission:
        model_object = transmission(atmosphereob)
        preob = preselector(model_object)

    elif params.fit_emission:
        model_object = emission(atmosphereob)
        preob = preselector(model_object)

    preob.run()
    logging.info('Molecule pre-selected using Marple module: %s' % preob.molselected)
    logging.info('Updating parameters object values')
    updated_params = preob.update_params()

    logging.info('Reinitialising data and atmosphereob objects')
    dataob = data(updated_params)
    #adding bulk composition to the atmosphere @todo move to better place
#     dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)
#     dataob.add_molecule('He', 4.0, 1.0e-9, 1.0000350, 0.15)
    
    atmosphereob = atmosphere(dataob)


#initialising transmission radiative transfer code instance
if params.fit_transmission:
    forwardmodelob = transmission(atmosphereob)

#initialising emission radiative transfer code instance
if params.fit_emission:
    forwardmodelob = emission(atmosphereob)

#initialising fitting object  @todo here we should actually have the transmission or emission objects as input. Then get rid of self.set_model()
fittingob = fitting(forwardmodelob)

#fit data
if params.downhill_run:
    # @todo should we run the downhill fit only on the first thread? Or maybe use different starting points
    fittingob.downhill_fit()    #simplex downhill fit

if params.mcmc_run and pymc_import:
    fittingob.mcmc_fit() # MCMC fit
    MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

if params.nest_run and multinest_import:
    fittingob.multinest_fit() # Nested sampling fit

#forcing slave processes to exit at this stage
if MPIimport and MPI.COMM_WORLD.Get_rank() != 0:
    #MPI.MPI_Finalize()
    exit()
    
    
#test area for cross correlation
T_mean, T_minmax, P = iterate_TP_profile(fittingob.NEST_FIT_mean[0], fittingob.NEST_FIT_std[0], atmosphereob.fit_index,atmosphereob.TP_profile)

T_sigma = T_minmax[:,1] - T_minmax[:,0]
Tmin = np.min(T_mean)
Tmax = np.max(T_mean)
nlayers = atmosphereob.nlayers

Trange = np.int(Tmax-Tmin + 100)
# print Trange

Tarray = np.linspace(Tmin-10,Tmax+10,Trange)
# Tarray = np.linspace(0,20,1000)

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / 2 * np.power(sig, 2.))


TPg_arr = np.zeros((nlayers,len(Tarray)))
Ent_arr = np.zeros((nlayers,nlayers))
Sig_arr = np.zeros((nlayers,nlayers))


# figure(100)
# 
# for i in range(nlayers):
#     meanoffset = T_mean[i] - Tmin + 100
#     TPg_arr[i,:] = gaussian(Tarray,T_mean[i],T_sigma[i])
# #     TPg_arr[i,:] = gaussian(Tarray,Tarray[500],T_sigma[i])
#     plot(Tarray,TPg_arr[i,:])
    


for i in range(nlayers):
    for j in range(nlayers):
#         Ent_arr[i,j] =  scipy.stats.entropy(pk=TPg_arr[i,:],qk=TPg_arr[j+i,:]) - scipy.stats.entropy(pk=TPg_arr[i,:])
#         Ent_arr[i,j] = - scipy.stats.entropy(pk=TPg_arr[i,:])
#         Ent_arr[i,j] = np.correlate(TPg_arr[i,:],TPg_arr[j,:])
#         Ent_arr[i,j] = np.abs(T_mean[i]/T_sigma[i]-T_mean[j]/T_sigma[j])
        Ent_arr[i,j] = np.abs(T_mean[i]-T_mean[j])
        Sig_arr[i,j] = np.sqrt(T_sigma[i]**2+T_sigma[j]**2)
#         Ent_arr[i,j] = T_sigma[i]/T_sigma[j+i]
    
    
covmatrix = np.zeros((nlayers,nlayers))
h=5.0
for i in xrange(nlayers):
    covmatrix[i,:] = np.exp(-1.0* np.abs(np.log(P[i]/P[:]))/h)      

Diff_arr = Ent_arr-Sig_arr

# figure(100)
# imshow(Diff_arr,origin='lower')
# colorbar()
# contour(Sig_arr,linewidths=2.0)
# 
# figure(101)
# imshow(Ent_arr,origin='lower')
# colorbar()
# contour(Sig_arr,linewidths=2.0)
# 
Diff_norm = ((Diff_arr-np.min(Diff_arr))/np.max(Ent_arr))
figure(102)
imshow(1.0-Diff_norm,origin='lower')
colorbar()
 
figure(103)
imshow(covmatrix,origin='lower')
colorbar()
contour(1.0-Diff_norm,linewidths=2.0)

figure(104)

subplot(141)
alpha = 1.0
imshow((1-alpha)*covmatrix+alpha*(1.0-Diff_norm),origin='lower')
contour(1.0-Diff_norm,linewidths=2.0)
title('alpha = 1.0')

subplot(142)
alpha = 0.7
imshow((1-alpha)*covmatrix+alpha*(1.0-Diff_norm),origin='lower')
contour(1.0-Diff_norm,linewidths=2.0)
title('alpha = 0.7')

subplot(143)
alpha = 0.4
imshow((1-alpha)*covmatrix+alpha*(1.0-Diff_norm),origin='lower')
contour(1.0-Diff_norm,linewidths=2.0)
title('alpha = 0.4')

subplot(144)
alpha = 0
imshow((1-alpha)*covmatrix+alpha*(1.0-Diff_norm),origin='lower')
# colorbar()
contour(1.0-Diff_norm,linewidths=2.0)
title('alpha = 0.0')


atmosphereob.set_TP_profile(profile='rodgers')

parble = [0.5,0.5]
for i in range(nlayers):
    parble.append(1.0)

# Th, Ph, ble = atmosphereob.TP_profile(parble,h=5.0)
# Tc,Pc,ble = atmosphereob.TP_profile(parble,covmatrix=(1.0-Diff_norm))
# 
# figure(103)
# plot(Th,Ph)
# plot(Tc,Pc,'k')
# yscale('log')
# xlabel('Temperature')
# ylabel('Pressure (Pa)')
# gca().invert_yaxis()
 
# figure(104)
# plot(T_sigma)

# print T_mean
# print T_sigma
show()
exit()
    
    

#initiating output instance with fitted data from fitting class
outputob = output(fittingob)

#plotting fits and data
logging.info('Plotting and saving results')

if params.verbose or params.out_save_plots:
    outputob.plot_all(save2pdf=params.out_save_plots)

# outputob.plot_spectrum()   #plotting data only
# outputob.plot_multinest()  #plotting multinest posteriors
# outputob.plot_mcmc()       #plotting mcmc posterios
# outputob.plot_fit()        #plotting model fits
#
outputob.save_model()       #saving models to ascii


#####################################################################
#launches external housekeeping script. E.g. useful to transfer data 
#from Scratch to home 
if params.clean_run:
    import subprocess
    #copies used parameter file to ./Output 
    if params.clean_save_used_params:
        subprocess.call('cp '+options.param_filename+' Output/',shell=True)

    subprocess.call('python '+params.clean_script,shell=True)


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
# if params.verbose:
#     show()
