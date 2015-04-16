#loading libraries
import pylab, sys, os, optparse, time, pickle
import numpy as np
from ConfigParser import SafeConfigParser
import math

sys.path.append('./classes')
import parameters
from parameters import *


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


RSOL  = 6.955e8         #stellar radius to m
RJUP  = 6.9911e7        #jupiter radius to m
MJUP  = 1.898e27        #jupiter mass to kg
REARTH= 6.371e3         #earth radius to m
AU    = 1.49e11         #semi-major axis (AU) to m
AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  help='Parameter filename')
options, remainder = parser.parse_args()

#initialising parameters object
params = parameters(options.param_filename)

# load first NEST_out to get fit params
for res in os.listdir(params.out_path):
    if os.path.isdir(os.path.join(params.out_path, res)):
        for error in os.listdir(os.path.join(params.out_path, res)):
            if os.path.isdir(os.path.join(params.out_path, res, error)):
                filename = os.path.join(params.out_path, res, error, 'NEST_out.db')
                if os.path.isfile(filename):
                    out = pickle.load(open(filename, 'rb'))
                    break
        break

# load a few things
fit_params = out[0]['fit_params']
X_names = params.planet_molec
X_in = params.planet_mixing
params_others = ['T', 'P0', 'Radius', 'mu', 'coupled_mu']

# print column labels
string = 'Resolution	Error (ppm)	'
for param in X_names:
    string += '%s	E(%s)	' % (param, param)
for param1 in X_names:
    for param2 in X_names:
        if param1 != param2:
            string += '%s/%s	E(%s/%s)	' % (param1, param2, param1, param2)
for param in out[0]['fit_params']:
    if not param in X_names:
        string += '%s	E(%s)	' % (param, param)

print string

# print input values
string = '		'
for idx, val in enumerate(X_names):
    string += '%s		' % (X_in[idx])
for idx1, val1 in enumerate(X_names):
    for idx2, val2 in enumerate(X_names):
        if val1 != val2:
            string += '%.4e		' % (X_in[idx1]/X_in[idx2])
for param in out[0]['fit_params']:
    if not param in X_names:
        if param == 'T':
            string += '%s		' % (params.planet_temp)
        if param == 'Radius':
            string += '%s		' % (params.planet_radius/RJUP)
        if param == 'P0':
            string += '%s		' % (params.tp_max_pres)
        if param == 'mu':
            string += '%s		' % (params.planet_mu/AMU)
print string

# print retrieval solutions
for res in os.listdir(params.out_path):
    if os.path.isdir(os.path.join(params.out_path, res)):
        for error in os.listdir(os.path.join(params.out_path, res)):
            if os.path.isdir(os.path.join(params.out_path, res, error)):
                filename = os.path.join(params.out_path, res, error, 'NEST_out.db')
                if os.path.isfile(filename):

                    out = pickle.load(open(filename, 'rb'))

                    for solution in out:

                        weights = solution['weights']

                        string = '%s	%s	' % (res[1:], error[1:])

                        # print molecules absolute abundances
                        for param in X_names:
                            if params.fit_X_log:
                                value = np.power(10, solution['fit_params'][param]['value'])
                                std = np.log(10) * value * solution['fit_params'][param]['std']
                            else:
                                value = solution['fit_params'][param]['value']
                                std = solution['fit_params'][param]['std']
                            string += '%.4e	%.4e	' % (value, std)

                        # print molecules abundance ratios
                        for param1 in X_names:
                            for param2 in X_names:
                                if param1 != param2:
                                    if params.fit_X_log:
                                        value1 = np.power(10, solution['fit_params'][param1]['value'])
                                        value2 = np.power(10, solution['fit_params'][param2]['value'])
                                        std1 = np.log(10) * value1 * solution['fit_params'][param1]['std']
                                        std2 = np.log(10) * value2 * solution['fit_params'][param2]['std']
                                    else:
                                        value1 = solution['fit_params'][param1]['value']
                                        value2 = solution['fit_params'][param2]['value']
                                        std1 = solution['fit_params'][param1]['std']
                                        std2 = solution['fit_params'][param2]['std']

                                    value = value1/value2
                                    std = value * np.sqrt((std1/value1)**2.+(std2/value2)**2.)

                                    string += '%.4e	%.4e	' % (value, std)

                        # print other params
                        for param in solution['fit_params']:
                            if not param in X_names:
                                if param == 'mu':
                                    value = solution['fit_params'][param]['value']/AMU
                                    std = solution['fit_params'][param]['std']/AMU
                                else:
                                    value = solution['fit_params'][param]['value']
                                    std = solution['fit_params'][param]['std']
                                string += '%.4e	%.4e	' % (value, std)

                        print string
