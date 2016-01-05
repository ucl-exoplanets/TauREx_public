'''

H2-H2 Collision Induced Absorption (400 - 7000 K) from A. Borysow

Create a pickled sigma array as a function of temperature for H2-H2 cia

Spectral range: 10 - 2000 wavenumbers at 10 wno resolution

Temperatures:  400, 500, 600, 700, 800, 900, 1000 - A. Borysow, A&A, v.390, p.779-782 (2002)
               1000, 2000, 3000, 4000, 5000, 6000, 7000 -  A. Borysow, U. G. Jorgensen and Y. Fu, , JQSRT, v.68, p.235-255 (2001)

'''

import numpy as np
import pickle

lowtemp_fname = 'Input/cia/data/final_CIA_400-1000.dat'
hightemp_fname = 'Input/cia/data/final_CIA_1000-7000.dat'

wngrid = np.arange(0, 50000, 10)

cia1 = np.loadtxt(lowtemp_fname)
cia2 = np.loadtxt(lowtemp_fname)

cia1_temp = [400, 500, 600, 700, 800, 900, 1000]
cia2_temp = [1000, 2000, 3000, 4000, 5000, 6000, 7000]

sigma_dict = {}
sigma_dict['t'] = np.asarray(cia1_temp + cia2_temp)

sigma_dict['xsecarr'] = np.zeros((len(cia1_temp+cia2_temp), len(wngrid)))
sigma_dict['wno'] = wngrid

for idx, val in enumerate(cia1_temp):
    sigma_dict['xsecarr'][idx] = np.interp(wngrid, cia1[:,0], cia1[:,idx+1])

for idx, val in enumerate(cia2_temp):
    sigma_dict['xsecarr'][len(cia1_temp) + idx] = np.interp(wngrid, cia2[:,0], cia1[:,idx+1])

sigma_dict['comments'] = ['400-1000K: A. Borysow, A&A, v.390, p.779-782 (2002)',
                          '1000-7000K: A. Borysow, U. G. Jorgensen and Y. Fu, , JQSRT, v.68, p.235-255 (2001) ']

pickle.dump(sigma_dict, open('Input/cia/H2H2.db', 'wb'))