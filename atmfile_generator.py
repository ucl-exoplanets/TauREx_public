#! /usr/bin/python 

# Small script that generates ATM files in the format: 
# Pressure (Pa), Temperature (K), Altitide (km) 
# #
# at the moment this is a stub but if needed i'll be made more user friendly

import numpy as np
import pylab as pl
import scipy.signal
from scipy import interpolate
import os
import scipy.constants as con


def find_nearest(arr, value):
    # find nearest value in array
    arr = np.array(arr)
    idx = (abs(arr-value)).argmin()
    return [arr[idx], idx]


def movingaverage(values,window):
    weigths = np.repeat(1.0, window)/window
    smas = np.convolve(values, weigths, 'valid')
#     smas2 = np.convolve(smas[::-1],weigths,'valid')
    return smas #smas2[::-1] # as a numpy array


N_LAYERS = 80 #number of atmospheric layers
N_SCALE  = 20.0   #number of scale heights 
MAX_P    = 1e6 #maximum pressure (Pa)

MOL_NUM  = 2   #number of molecules 

smooth_window = 5 #smoothing window size as percent of total data

T_surf = 1500.0  #surface temperature (K)
mu     = 2.3   #mean molecular weight (atomic units)
g      = 7.1  #surface gravity (m/s^2)
        

# print 'boltzmann ', con.k
# print 'T_surf ', T_surf
# print 'mmw ', mu * 1.660538921e-27
# print 'g ',g

scaleheight = (con.k * T_surf) / ((mu * 1.660538921e-27) * g)

# print scaleheight

max_z    = N_SCALE * scaleheight

PTA_arr = np.zeros((N_LAYERS,3+MOL_NUM))
PTA_arr[:,2] = np.linspace(0,max_z,num=N_LAYERS)
PTA_arr[:,0] = MAX_P * np.exp(-PTA_arr[:,2]/scaleheight)
PTA_arr[:,1] = T_surf

MIN_P = PTA_arr[-1,0]
print 'min_p ', PTA_arr[-1,0]

print 'scaleheight ', scaleheight
print 'N_SCALE ', N_SCALE
print 'max_z ',max_z, 'end z ', PTA_arr[-1,2]


#setting TP nodes
# Znodes = [0.0, 4.0e6,5e6,max_z]
# Pnodes = [MAX_P,1e5, 100.0,MIN_P]
Pnodes = [MAX_P, 1000.0,MIN_P]
# Tnodes = [T_surf,1200.0,1200.0]
# Tnodes = [T_surf,T_surf,1200.0,1200.0]
Tnodes = [T_surf,1200.0,1200.0]
# Tnodes = [T_surf,T_surf,T_surf, T_surf]

print 'Pnodes ',Pnodes

#creating linear T-P profile
# int_func1 = interpolate.interp1d(Pnodes[::-1],Tnodes[::-1],kind='quadratic')
# TP = int_func1((PTA_arr[::-1])[:,0])


TP = np.interp((np.log(PTA_arr[::-1])[:,0]), np.log(Pnodes[::-1]), Tnodes[::-1])

# print TP
#smoothing T-P profile
wsize = N_LAYERS*(smooth_window/50.0)
if (wsize %2 == 0):
    wsize += 1 
TP_smooth = movingaverage(TP,wsize)
border = np.int((len(TP) - len(TP_smooth))/2)


PTA_arr[:,1] = TP[::-1]
PTA_arr[border:-border,1] = TP_smooth[::-1]



#adding mixing ratios for molecules 
# X = [6.0e-5,1.0e-5,5.0e-6]
# X = [1.0e-7,1.0e-7,1.0e-7]
X = [1e-6,2e-6]

for i in range(len(X)):
    PTA_arr[:,2+i+1] = X[i]
    


np.savetxt('TP_profile.txt', PTA_arr)


pl.figure(1)
# pl.plot(TP, PTA_arr[:,0],'k')
pl.plot(PTA_arr[:,1],PTA_arr[:,0],'r')
pl.plot(Tnodes,Pnodes, 'x')
pl.xlim(np.min(Tnodes)-np.min(Tnodes)*0.1,np.max(Tnodes)+np.max(Tnodes)*0.1)
pl.yscale('log')
pl.xlabel('Temperature')
pl.ylabel('Pressure (Pa)')
pl.gca().invert_yaxis()

# pl.figure(2)
# pl.plot(PTA_arr[:,0])
# pl.title('pressure')
# pl.figure(3)
# pl.plot(PTA_arr[:,1],'r')
# pl.title('temperature')
# pl.figure(4)
# pl.plot(PTA_arr[:,2],'g')
# pl.title('z')


pl.show()
# pl.draw()