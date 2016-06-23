'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Convert k tables to a new grid of pressures and temperatures

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


parser = argparse.ArgumentParser()

parser.add_argument('--input_file', '--input_filename',
                  dest='input_file',
                  default=False,
)
parser.add_argument('--output_file', '--output_filename',
                  dest='output_file',
                  default=False,
)
parser.add_argument('-p', '--pressure_list',
                  dest='pressure_list',
                  default=False,
                  )
parser.add_argument('-t', '--temperature_list',
                  dest='temperature_list',
                  default=False,
                  )

options = parser.parse_args()

if not options.input_file:
    print 'You need to specify an input file or folder'
    exit()
if not options.output_file:
    print 'You need to specify an output file'
    exit()

ktable = pickle.load(open(options.input_file))


if not options.pressure_list:
   pressures = ktable['p']
else:
   l = options.pressure_list.split(',')
   pressures = np.asarray([float(m) for m in l])

if not options.temperature_list:
   temperatures = ktable['t']
else:
    l = options.temperature_list.split(',')
    temperatures = np.asarray([float(m) for m in l])

pressures = np.sort(pressures)
temperatures = np.sort(temperatures)

kcoeff = ktable['kcoeff']

if options.temperature_list:
    kcoeff_tmp_1 = np.zeros((len(ktable['p']), len(temperatures), len(ktable['wno']), ktable['ngauss']))
    for pressure_idx, pressure_val in enumerate(ktable['p']):
        for temperature_idx, temperature_val in enumerate(temperatures):
            for wno_idx, wno_val in enumerate(ktable['wno']):
                for kcoeff_idx in range(ktable['ngauss']):
                    kcoeff_tmp_1[pressure_idx, temperature_idx, wno_idx, kcoeff_idx] = \
                        np.interp(np.log10(temperature_val), np.log10(ktable['t']), kcoeff[pressure_idx, :, wno_idx, kcoeff_idx])

else:
    kcoeff_tmp_1 = kcoeff

if options.pressure_list:
    kcoeff_tmp_2 = np.zeros((len(pressures), len(temperatures), len(ktable['wno']), ktable['ngauss']))
    for pressure_idx, pressure_val in enumerate(pressures):
        for temperature_idx, temperature_val in enumerate(temperatures):
            for wno_idx, wno_val in enumerate(ktable['wno']):
                for kcoeff_idx in range(ktable['ngauss']):
                    kcoeff_tmp_2[pressure_idx, temperature_idx, wno_idx, kcoeff_idx] = \
                        np.interp(np.log10(pressure_val), np.log10(ktable['p']), kcoeff_tmp_1[:, temperature_idx, wno_idx, kcoeff_idx])

else:
    kcoeff_tmp_2 = kcoeff_tmp_1

ktable_out = {}
ktable_out = {}
ktable_out['wno'] = ktable['wno']
ktable_out['resolution'] = ktable['resolution']
ktable_out['weights'] = ktable['weights']
ktable_out['ngauss'] = ktable['ngauss']
ktable_out['kcoeff'] = kcoeff_tmp_2
ktable_out['name'] = ktable['name']
ktable_out['p'] = pressures
ktable_out['t'] = temperatures

pickle.dump(ktable_out, open(options.output_file, 'wb'), protocol=2)