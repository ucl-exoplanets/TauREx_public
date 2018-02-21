#! /usr/bin/python 

#small script that shits out the venot file format equivalent for the lightening project

import numpy as np
import pylab as pl
import pyfits as pf

import glob, os, sys


AMU = 1.660538921e-27
KBOLTZ = 1.380648813e-23
G = 6.67384e-11
RSOL = 6.955e8
RJUP = 6.9911e7
#RJUP = 7.1492e7 # Jo's radius
MJUP = 1.898e27
AU = 1.49e11




DIR = '/Users/ingowaldmann/Dropbox/UCLlocal/REPOS/taurex/Input/lightening'
# FILENAME = 'Earth-Today-Lightning-Full.dat'
FILENAME = 'Modern-Earth-noLightning-Full.dat'

#EARTH
planet_mass = 5.97237e24 #kg
planet_radius = 6.371e6 #m
planet_mu   =  28.97 * AMU#kg

data = np.loadtxt(os.path.join(DIR,FILENAME),skiprows=1)
[nlayers,ncols] = np.shape(data)

fheader = open(os.path.join(DIR,FILENAME),'r')
header = fheader.readlines()
c=0
for line in header:
    head = line
    break
fheader.close()

#rebuilding header line 
newhead = 'alt(km)  '+head[:2]+'m'+head[2:]
newhead_small = head[62:]

print head.split()


molnames = ['C_1D','H','N','O','O_1D','O_1S', 'CO', 'H2','HO','N2', 'NO', 'O2', 'O2_D', 'O3', 'CH4', 
            'CO2', 'H2O', 'HO2', 'N2O', 'NO2', 'H2O2', 'HNO3', 'CH2O2', 'HCOOH', 'CH3ONO', 'e-', 
            'H+', 'O+', 'NO+','O2+', 'C','HN','CNC','H2N','H3N','C+','C-','N+','O-','CO+','HO+','N2+','CHO+',
            'CH3','CHO','HCN','HNO','NO3','C2H2','C2H6','CH2O','HNO2','N2O3','CH3O2','CH3OH','CH4O2','H3O+']
molweights = [14,1,14,16,18,18,28,2,17,28,30,32,34,48,16,44,18,33,44,46,34,63,46,46,61,0,1,16,30,32,12,15,38,16,17,12,12,14,16,28,17,28,29,
              15,29,27,31,62,26,28,30,48,76,47,32,52,19]


badwords = ['p(bar)' ,'T(K)' , 'NH(cm-3)' , 'Kzz(cm2s-1)' , 'Hz(cm)', 'zeta(s-1)']

mollist = []
for mol in head.split():
    if mol in molnames:
        mollist.append(molweights[molnames.index(mol)])
    
    elif mol not in badwords:
        mollist.append(2.3)
        print 'FILLED: ',mol
    else:
        print 'OUT: ',mol
#         mollist.append(2.3)

moleweigthstr = ' '.join(str(e) for e in mollist)

#create ranking of most important molecules according to abundance
molabundance =[]
mnamelist =[]
c=0
for mol in head.split():
    mnamelist.append(mol)
    molabundance.append(np.max(data[:,c]))
    c+=1

mnamelist = np.asarray(mnamelist[6:])
molabundance = np.asarray(molabundance[6:])

midx = np.argsort(molabundance)

print midx[::-1]
print mnamelist[midx][::-1]
print molabundance[midx][::-1]


pressure_profile_levels = data[:,0] * 1000.0 #converting bar to mbar 
temperature_profile = data[:,1]


H = np.zeros(nlayers)
g = np.zeros(nlayers)
z = np.zeros((nlayers,1))

g[0] = (G * planet_mass) / (planet_radius**2) # surface gravity (0th layer)
H[0] = (KBOLTZ*temperature_profile[0])/(planet_mu*g[0]) # scaleheight at the surface (0th layer)


for i in xrange(1, nlayers):
    deltaz = (-1.)*H[i-1]*np.log(pressure_profile_levels[i]/pressure_profile_levels[i-1])
    z[i] = z[i-1] + deltaz # altitude at the i-th layer
    with np.errstate(over='ignore'):
        g[i] = (G * planet_mass) / ((planet_radius + z[i])**2) # gravity at the i-th layer
    with np.errstate(divide='ignore'):
        H[i] = (KBOLTZ*temperature_profile[i])/(planet_mu*g[i])


z /=1e3 #converting m to km

OUT = np.hstack((z,data))
OUT2 = OUT[:,:3]

[s1,s2] = np.shape(data[:,6:])
OUT3 = np.zeros((s1,s2+2))
OUT3[:,0] = z[:,0]
OUT3[:,1] = data[:,0]
OUT3[:,2:] = data[:,6:]



with open(FILENAME[:-4]+'_conv.dat','wb') as outfile: 
    outfile.write(newhead)
    outfile.write(moleweigthstr+'\n')
    np.savetxt(outfile, OUT)
    
with open(FILENAME[:-4]+'_mixing.dat','wb') as outfile: 
    outfile.write(newhead_small)
    outfile.write(moleweigthstr+'\n')
    np.savetxt(outfile, OUT3)

np.savetxt(FILENAME[:-4]+'_tp.dat',OUT2)

pl.figure(1)
pl.plot(np.log(molabundance[midx][::-1]),linewidth=3.0)
pl.gca().xaxis.set_ticks(np.arange(0, len(molabundance), 1.0))
pl.gca().set_xticklabels(mnamelist[midx][::-1])
pl.ylabel('log(mixing ratio)')

pl.show()