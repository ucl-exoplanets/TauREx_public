'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Compute CIA cross sections

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''
import numpy as np
import pickle
from time import gmtime, strftime

'''

from HITRAN

Wavenumber range: 20 - 10000 in steps of 1 wn
T range: 200 - 3000 K in steps of 25K

Reference: C. Richard, I.E. Gordon, L.S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, et al, JQSRT 113, 1276-1285 (2012

Cross sections are in units of cm^-1 mol^-2. Need to convert to (m^5 mol^-2):
 convert from cm^5 to m^5: factor 1e-10

'''


# H2-H2

fname = 'Input/cia/HITRAN/data/H2-H2_2011.cia'

lines = [line.rstrip('\n') for line in open(fname)]

wngrid = np.arange(20, 10001, 1)
sigma_dict = {}
sigma_dict['wno'] = wngrid
sigma_dict['t'] = np.asarray([float(line_val.split()[4]) for line_idx, line_val in enumerate(lines) if line_idx == 0 or line_idx%9982 == 0.])

sigma_dict['xsecarr'] = np.zeros((len(sigma_dict['t']), len(wngrid)))

xsec = []
k = 0
for line_idx, line_val in enumerate(lines):
    if line_idx == 0:
        continue
    elif line_idx%9982 == 0.:
        sigma_dict['xsecarr'][k] = np.asarray(xsec)
        xsec = []
        k += 1
    else:
        xsec.append(float(line_val.split()[1])*1e-10) # convert from cm^5 to m^5
sigma_dict['comments'] = ['H2-H2 collision induced absorption',
                          '200-3000K in steps of 25K',
                          'C. Richard, I.E. Gordon, L.S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, et al, JQSRT 113, 1276-1285 (2012).',
                          'Created with TauREx at GMT %s' % strftime("%Y-%m-%d %H:%M:%S", gmtime())]
pickle.dump(sigma_dict, open('Input/cia/HITRAN/H2-H2.db', 'wb'))

# H2-He

fname = 'Input/cia/HITRAN/data/H2-He_2011a.cia'

lines = [line.rstrip('\n') for line in open(fname)]

wngrid = np.arange(20, 20001, 1)
sigma_dict = {}
sigma_dict['wno'] = wngrid
sigma_dict['t'] = np.asarray([float(line_val.split()[4]) for line_idx, line_val in enumerate(lines) if line_idx == 0 or line_idx%19982 == 0.])
sigma_dict['t'] = sigma_dict['t'][sigma_dict['t']<=3000] # exclude temperatures higher than 3000 K

sigma_dict['xsecarr'] = np.zeros((len(sigma_dict['t']), len(wngrid)))

xsec = []
k = 0
for line_idx, line_val in enumerate(lines):
    if line_idx == 0:
        continue
    elif line_idx%19982 == 0.:
        sigma_dict['xsecarr'][k] = np.asarray(xsec)
        xsec = []
        k += 1
        if float(line_val.split()[4]) > 3000: # exit loop for temperatures higher than 3000 K
            break
    else:
        xsec.append(float(line_val.split()[1])*1e-10) # convert from cm^5 to m^5
sigma_dict['comments'] = ['H2-He collision induced absorption',
                          '200-3000K in steps of 25K',
                          'C. Richard, I.E. Gordon, L.S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, et al, JQSRT 113, 1276-1285 (2012).',
                          'Created with TauREx at GMT %s' % strftime("%Y-%m-%d %H:%M:%S", gmtime())]

pickle.dump(sigma_dict, open('Input/cia/HITRAN/H2-HE.db', 'wb'))

exit()

'''

from A. Borysow

Compute CIA cross sections. Assume the same temperature grid for all pairs.
Note that for H2-He there are no cross sections for T < 1000K, assume zero.

The temperature grid (in K) is 400, 500, 600, 700, 800, 900, 1000, 1000, 2000, 3000, 4000, 5000, 6000, 7000

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
H2-H2 Collision Induced Absorption (400 - 7000 K) from A. Borysow

Create a pickled sigma array as a function of temperature for H2-H2 cia

Spectral range: 10 - 2000 wavenumbers at 10 wno resolution

Temperatures:  400, 500, 600, 700, 800, 900, 1000 - A. Borysow, A&A, v.390, p.779-782 (2002)
               1000, 2000, 3000, 4000, 5000, 6000, 7000 -  A. Borysow, U. G. Jorgensen and Y. Fu, , JQSRT, v.68, p.235-255 (2001)


Cross sections are in units of cm^-1 amagat^-2. Need to convert to (m^5 mol^-2):
 convert from (cm^-1 amagat^-2) to (cm^5 mol^-2): factor = 1.0 / np.power((AMAGAT*1.0e-6),2) = 1.385294e-39
 convert from cm^5 to m^5: factor 1e-10

'''

lowtemp_fname = 'Input/cia/Borysow/data/final_CIA_400-1000.dat'
hightemp_fname = 'Input/cia/Borysow/data/final_CIA_1000-7000.dat'

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
    sigma_dict['xsecarr'][len(cia1_temp) + idx] = np.interp(wngrid, cia2[:,0], cia2[:,idx+2])


sigma_dict['comments'] = ['400-1000K: A. Borysow, A&A, v.390, p.779-782 (2002)',
                          '2000-7000K: A. Borysow, U. G. Jorgensen and Y. Fu, , JQSRT, v.68, p.235-255 (2001)',
                          'Created with TauREx at GMT %s' % strftime("%Y-%m-%d %H:%M:%S", gmtime())]

pickle.dump(sigma_dict, open('Input/cia/Borysow/H2-H2.db', 'wb'))


# h2-he pair

#lowtemp_fname = 'Input/cia/data/?'
hightemp_fname = 'Input/cia/data/Borysow/h2he_CIA_1000-7000.dat'

wngrid = np.arange(0, 50000, 10)

#cia1 = np.loadtxt(lowtemp_fname)
cia2 = np.loadtxt(hightemp_fname)

cia1_temp = [400, 500, 600, 700, 800, 900]
cia2_temp = [1000, 2000, 3000, 4000, 5000, 6000, 7000]

sigma_dict = {}
sigma_dict['t'] = np.asarray(cia1_temp + cia2_temp)

sigma_dict['xsecarr'] = np.zeros((len(cia1_temp+cia2_temp), len(wngrid)))
sigma_dict['wno'] = wngrid

# for idx, val in enumerate(cia1_temp):
#     sigma_dict['xsecarr'][idx] = np.interp(wngrid, cia1[:,0], cia1[:,idx+1])

for idx, val in enumerate(cia2_temp):
    sigma_dict['xsecarr'][len(cia1_temp) + idx] = np.interp(wngrid, cia2[:,0], cia2[:,idx+1])*1.385294e-49

sigma_dict['comments'] = ['400-900K: N/A, assume zero',
                          '1000-7000K: A. Borysow, U. G. Jorgensen and C. Zheng, A&A, vol. 324, p.185-195 (1997)',
                          'Created with TauREx at GMT %s' % strftime("%Y-%m-%d %H:%M:%S", gmtime())]

pickle.dump(sigma_dict, open('Input/cia/Borysow/H2-HE.db', 'wb'))