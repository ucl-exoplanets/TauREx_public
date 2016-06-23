'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Collect k tables

    Developers: Marco Rocchetto (University College London)

'''

import numpy as np
import cPickle as pickle
import os
from ConfigParser import SafeConfigParser
import argparse
import time
import glob
import sys

#loading parameter file parser
parser = argparse.ArgumentParser()

parser.add_argument('--input_files',
                      dest='input_files',
                      default=False,
                    )
parser.add_argument('--output_file',  # can be pickle or ascii
                      dest='output_file',
                      default=False
                    )

options = parser.parse_args()

if not options.input_files:
    print 'You need to specify an input file or folder'
    exit()
if not options.output_file:
    print 'You need to specify an output file'
    exit()


files = glob.glob(os.path.join(options.input_files))
ktables = []

for filenm in files:
    ktable = pickle.load(open(filenm))

    ktables.append(ktable)

wno  = ktables[0]['wno']
resolution  = ktables[0]['resolution']
weights  = ktables[0]['weights']
ngauss  = ktables[0]['ngauss']
name  = ktables[0]['name']
wlrange  = ktables[0]['wlrange']
wnrange  = ktables[0]['wnrange']
method  = ktables[0]['method']

pressures = []
temperatures = []

for idx, ktable in enumerate(ktables):
    if len(ktable['wno']) != len(wno) or ktable['resolution'] != resolution or \
    len(ktable['weights']) != len(weights) or ktable['ngauss'] != ngauss or ktable['name'] != name or \
     ktable['wlrange'] != wlrange or ktable['wnrange'] != wnrange or ktable['method'] != method:
        print 'Something is not right with %s' % files[idx]
        exit()

    pressures.append(ktable['p'])
    temperatures.append(ktable['t'])

pressures = np.sort(np.unique(pressures))
temperatures = np.sort(np.unique(temperatures))

kcoeff = np.zeros((len(pressures), len(temperatures), len(wno), ngauss))

for idx, ktable in enumerate(ktables):

    p_idx = np.where(pressures==ktable['p'])[0][0]
    t_idx = np.where(temperatures==ktable['t'])[0][0]

    kcoeff[p_idx, t_idx, :, :] = ktable['kcoeff']

ktable_out = {}
ktable_out['wno'] = wno
ktable_out['resolution'] = resolution
ktable_out['weights'] = weights
ktable_out['ngauss'] = ngauss
ktable_out['kcoeff'] = kcoeff
ktable_out['name'] = name
ktable_out['t'] = temperatures
ktable_out['p'] = pressures
ktable_out['wlrange'] = wlrange
ktable_out['wnrange'] = wnrange
ktable_out['method'] = method

pickle.dump(ktable_out, open(options.output_file, 'wb'), protocol=2)