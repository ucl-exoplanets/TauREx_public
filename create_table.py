'''

Creates a table (CSV format) with the retrieved parameters for a given run of create_grid.py. Reads in grid.db

Format: param_name
        input_value
        resolution
        snr
        mode (default = 0)
        [retrieved_param, retrieved_param_err, relative_error]  x modes

        where modes is the number of modes found by multinest

'''

import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib import rc
import optparse
import pymultinest
import os
import sys
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5 #set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 8})


sys.path.append('./classes')
sys.path.append('./library')
import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

parser = optparse.OptionParser()
parser.add_option('-i', '--input',
                  dest='input',
                  default=False,
)
parser.add_option('-o', '--output',
                  dest='output',
                  default='grid.csv',
)
parser.add_option('-n', '--name',
                  dest='name',
                  default=False,
)

options, remainder = parser.parse_args()

if not options.input or not options.output:
    print parser.print_help()
    exit()

# load database
db = pickle.load(open(options.input, 'rb'))

plt.savefig(os.path.join(os.path.split(options.input)[0], 'contour_plots.pdf'))


# Build table
grid_filename = options.output

f = open(grid_filename,'w')
for res in db['resolutions']: # loop over molecules
    for snr in db['snrs']:
        fit = db['model_fitting'][res][snr]

        line = 'Temperature, %i, %i, %i, 0' %  (db['temperature_input'],res, snr)
        modes = fit['modes']
        for mode in range(modes):
            line += ', %.2e, %.2e, %.2e' % (fit['T'][mode], fit['T_std'][mode],
                                            fit['T_std'][mode]/fit['T'][mode])
        line += '\n'
        f.write(line)

for mol in range(len(db['molecules_input'])): # loop for each molecule
    for res in db['resolutions']:
        for snr in db['snrs']:
            fit = db['model_fitting'][res][snr]
            line = '%s, %.2e, %i, %i, 0' % (db['molecules_input'][mol], db['mixing_input'][mol], res, snr)
            modes = fit['modes']
            for mode in range(modes):
                line += ', %.2e, %.2e, %.2e' % (fit['X'][mode][mol][0], fit['X_std'][mode][mol],
                                                fit['X_std'][mode][mol]/fit['X'][mode][mol][0])
            line += '\n'
            f.write(line)