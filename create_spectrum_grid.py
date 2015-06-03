#! /usr/bin/python -W ignore

###########################################################
#    
#       
###########################################################

#loading libraries     
import numpy, pylab, sys, os, optparse, time,logging,itertools,gzip
from numpy import * #nummerical array library 
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser

from multiprocessing import Pool
# from itertools import product
from functools import partial

import cPickle as pickle
from library.library_emission import fit_brightness_temp

#checking for multinest library
# global multinest_import 
# try: 
#     with os.devnull as sys.stdout: 
#         import pymultinest
#         multinest_import = True
# except:
#     multinest_import = False

#loading classes
sys.path.append('./classes')
sys.path.append('./library')

import parameters,emission,transmission,fitting,atmosphere,data
from parameters import *
from emission import *
from transmission import *
from fitting import *
from atmosphere import *
from data import *

#loading libraries
import library_emission, library_transmission, library_general, library_plotting
from library_emission import *
from library_transmission import *
from library_general import *
from library_plotting import *

import create_spectrum 
from create_spectrum import *


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
parser.add_option('-s', '--save',
                  dest="save_name",
                  default="spectrum_grid"
)
parser.add_option('-n', '--n_processes',
                  dest="Nproc",
                  default=4
)

options, remainder = parser.parse_args()

#####################################################################

#constants
RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg

#parameters
N_processes =  int(options.Nproc)  #number of cores on which to run this
N_xsteps    = 7         #number of steps in X grid
N_tpsteps   = 7         #number of steps in TP grid


#settings some arrays, dics and lists
tpparlist = [] #list of tp parameter iterations
xparlist  = [] #list of X parameter iterations
bulk      = [] #list of bulk planetary properties: Rp, Mp, Ts, Rs
Pardic    = {} #dictionary that will hold all inputs and outputs



#X bounds for molecules
X_bounds = [(1e-8,1e-3),(1e-8,1e-3),(1e-8,1e-3),(1e-8,1e-3)]
#Bounds for '2-point' TP profile
TP_bounds = [(1000.0,2500.0), (100.0,700.0),(1e2,1e3)]


#setting bulk parameter grid: Rp, Mp, Ts, Rs
#doing this manually now since I only want a few
bulk.append([1.138,1.138,4800.0,0.8])  #hd 189733b type
bulk.append([0.24,0.020,3026.0,0.216]) #gj 1214b type
bulk.append([1.736,1.404,6300.0,1.6])  #wasp 12b type
bulk.append([0.422,0.082,4750.0,0.81]) #hat-p-11b type
bulk.append([1.38,0.69,6000.0,1.203])  #hd 209458b type



#filling TP and X parameter iteration arrays
for bound in X_bounds:
    xparlist.append(np.linspace(bound[0], bound[1], N_xsteps))    
xparlist = np.asarray(xparlist)

for bound in TP_bounds:
    tpparlist.append(np.linspace(bound[0], bound[1], N_tpsteps))

#calculating possible TP parameter permutations
tp_variations = np.array([x for x in itertools.product(*tpparlist)])
N_tpvar = np.shape(tp_variations)[0]    #number of permutations 


#setting more parameters
N_planets = len(bulk)                   #number of planets considered
N_gas     = len(X_bounds)               #number of gases 
N_tpvar   = np.shape(tp_variations)[0]  #number of permutations 


#timing 
start_time = time.clock()

#generating Pardic: [planet index][spectrum index][content]
count = 0
for p in range(N_planets):
    Pardic[p] = {}
    count = 0
    for i in range(N_gas):
        for j in range(N_xsteps):
            for k in range(N_tpvar):
                Pardic[p][count] = {}
                tmp = np.zeros(N_gas)
                tmp[i] = xparlist[i,j]
                Pardic[p][count]['X']        = tmp
                Pardic[p][count]['TP']       = tp_variations[k]
                Pardic[p][count]['T_planet'] = tp_variations[k][0]
                Pardic[p][count]['R_planet'] = bulk[p][0] * RJUP
                Pardic[p][count]['M_planet'] = bulk[p][1] * MJUP
                Pardic[p][count]['T_star']   = bulk[p][2] 
                Pardic[p][count]['R_star']   = bulk[p][3] * RSOL
                count += 1


#getting number of iterations
N_iter       = len(Pardic[0].keys())  #iterations per planet
N_iter_total = N_iter*N_planets       #total iterations

           
           
#loading create_spectrum            
#defining options for create_spectrum
options.param_filename = 'Parfiles/taurex_emission_wasp76_grid.par' #needs to be a valid parameter file but not important which
options.error = 0.0
options.bin = "resolution"
options.resolution = 100
options.noise = 0.0
#loading object
createob = create_spectrum(options)
logging.getLogger('').disabled = True #disables any further Tau-REx print statements 
params = createob.params              #making a copy of params object 


#function called by pool.map() to compute spectrum, normalise spectrum and return 
def populate_model_array(p,idx): 
    
    #setting bulk parameters
    createob.fmob.atmosphere.T[:] = Pardic[p][idx]['T_planet'] #set bulk planetary temperature from 10bar temperature 
    #setting abundances
    createob.set_mixing_ratios(Pardic[p][idx]['X'])
    #setting TP profile 
    createob.generate_tp_profile_2(Pardic[p][idx]['TP']) #set TP-profile
    #generating spectrum
    
    model = createob.generate_spectrum()    #generating spectrum
    model = model[~np.isnan(model).any(1)]
    tp_profile = createob.fmob.atmosphere.T
    model_wave = model[:,0] 
    
#    print np.min(tp_profile),np.max(tp_profile)

    bb_planet_l = black_body(model[:,0],np.min(tp_profile))#*0.88
    bb_planet_h = black_body(model[:,0],np.max(tp_profile))
    bb_star     = black_body(model[:,0],Pardic[p][idx]['T_star'])
#     bb_star     = black_body(model[:,0],4800)
    
#     bb_star *= (Pardic[p][idx]['R_star']**2*4.0*np.pi)
    
    F_star = createob.fmob.F_star    
    FpFs_l = bb_planet_l/bb_star *(Pardic[p][idx]['R_planet']/Pardic[p][idx]['R_star'])**2
    FpFs_h = bb_planet_h/bb_star *(Pardic[p][idx]['R_planet']/Pardic[p][idx]['R_star'])**2
    
#    figure()
#    plot(model_wave,model[:,1])
#    plot(model_wave,FpFs_l,'g')
#    plot(model_wave,FpFs_h,'r')  

    #normalising with blackbodies
    norm1 = model[:,1] - FpFs_l
    norm1 /= (FpFs_h-FpFs_l)
    
    #normalising first and second moment
    norm2 = norm1 - np.mean(norm1)
    norm2 /= np.std(norm2)
    
#     #fitting for brightness temperature 
#     tempgrid = np.zeros_like(model[:,1])
#     
#     for i in xrange(len(tempgrid)):
#         tempgrid[i] = fit_brightness_temp(model[i,0],model[i,1])
    
#     figure()
#     plot(bb_planet_l,'g')
#     plot(bb_planet_h,'r')

#     figure()
#     plot(createob.fmob.specgrid,F_star,'r')
#     plot(model_wave,bb_star,'k')
#     
#     figure()
#     plot(model_wave,norm2)
    
#     figure()
#     plot(model_wave,tempgrid)
    
    return model_wave, model[:,1], norm2 #returns: wavelength grid, spectrum, normalised spectrum


#setting up additional output arrays
testwave,test,norm_test = populate_model_array(0,0) #just getting array sizes
modelarray = np.zeros((len(testwave),N_iter_total)) #non dictionary output array containing the normalised model

# show()
# exit()


#do not delete following line. useful some day. some fine day. 
# iteridx = [(i,j) for i, j in itertools.product([p for p in range(N_planets)],[i for i in range(10)])]


#populating dictionary and modelarray
for p in range(N_planets):
    print 'Planet No.: ',p
    #whenever new planet is processed, reset createob with new bulk parameters
    params.planet_radius = Pardic[p][0]['R_planet']
    params.planet_mass   = Pardic[p][0]['M_planet']
    params.star_temp     = Pardic[p][0]['T_star']
    params.star_radius   = Pardic[p][0]['R_star']
    createob.reset(options, params)

    
    func = partial(populate_model_array, p) #wrapper to accept extra inputs to pool.map()
    pool = Pool(processes=N_processes)      #setting number of cores on which to run
    pool_result = pool.map(func,[i for i in range(N_iter)]) #runnning the stuff
    pool.close()                            #closing the pool 
    pool.join()                             #closing the pool. Needs to be done otherwise createob.reset has no effect.
    
    #collecting results and sorting them into dictionary and modelarray
    for idx, result in enumerate(pool_result):
        wave,model,norm = result 
        Pardic[p][idx]['wave']        = wave
        Pardic[p][idx]['model']       = model
        Pardic[p][idx]['norm']        = norm
        modelarray[:,idx+(N_iter*p)]  = norm



#saving dictionary
with gzip.open(options.save_name+'.pkl.gz', 'wb') as handle:
    pickle.dump(Pardic, handle)
np.savetxt(options.save_name+'.gz',modelarray)

end_time = time.clock()

print 'Total number of spectra: ', N_iter_total
print 'Total runtime: %.2fm' % ((end_time - start_time)/ 60.)

# figure()
# for i in range(200):
#     plot(modelarray[:,i])
# # plot(modelarray[:,15])

# figure()
# for i in range(100):
#     plot(Pardic[0][i]['model'])
# 
# figure()
# for i in range(100):
#     plot(Pardic[1][i]['model'])
# show()

