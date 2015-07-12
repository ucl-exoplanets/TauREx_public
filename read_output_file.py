#loading libraries
import pylab, sys, os, optparse, time, pickle
import numpy as np
from ConfigParser import SafeConfigParser
import math

AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-i', '--input',
                  dest="db_file",
)

options, remainder = parser.parse_args()
out = pickle.load(open(options.db_file))

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


params_clr = ['14N-1H3_CLR', '1H2-16O_CLR', '12C-16O2_CLR', '12C-1H4_CLR', '1H2-32S_CLR', '14N-1H3_CLR', 'H2_CLR', 'N2_CLR']
params_gas = ['14N-1H3', '1H2-16O', '12C-16O2', '12C-1H4', '1H2-32S', '14N-1H3', 'H2', 'N2']
params_others = ['T', 'P0', 'Radius', 'coupled_mu']

for solution in out:
    string = ''

    # for param in params_gas:
    #     value = solution['fit_params'][param]['value']
    #     std = weighted_avg_and_std(solution['fit_params'][param]['trace'], solution['weights'])
    #     string += '%.4e	%.4e	' % (value, std)

    for param in params_gas:
        value = np.power(10, solution['fit_params'][param]['value'])
        std = weighted_avg_and_std(np.power(10,solution['fit_params'][param]['trace']), solution['weights'])[1]
        string += '%s	%.4e	%.4e	' % (param, value, std)

    # # H20/CH4
    # value = np.power(10, solution['fit_params']['1H2-16O']['value']-solution['fit_params']['12C-1H4']['value'])
    # std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']-solution['fit_params']['12C-1H4']['trace']),
    #                            solution['weights'])[1]
    # string += '%.4e	%.4e	' % (value, std)
    #
    # # H2O/NH3
    # value = np.power(10, solution['fit_params']['1H2-16O']['value']-solution['fit_params']['14N-1H3']['value'])
    # std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']-solution['fit_params']['14N-1H3']['trace']),
    #                            solution['weights'])[1]
    # string += '%.4e	%.4e	' % (value, std)
    #
    # # CH4/NH3
    # value = np.power(10, solution['fit_params']['12C-1H4']['value']-solution['fit_params']['14N-1H3']['value'])
    # std = weighted_avg_and_std(np.power(10, solution['fit_params']['12C-1H4']['trace']-solution['fit_params']['14N-1H3']['trace']),
    #                            solution['weights'])[1]
    # string += '%.4e	%.4e	' % (value, std)

    # #He+H2+N2
    # value = np.power(10, solution['fit_params']['He']['value'])+ np.power(10, solution['fit_params']['H2']['value']) + \
    #         np.power(10, solution['fit_params']['N2']['value'])
    # std = weighted_avg_and_std(np.power(10, solution['fit_params']['He']['trace']) +
    #                            np.power(10, solution['fit_params']['H2']['trace']) +
    #                            np.power(10, solution['fit_params']['N2']['trace']),
    #                            solution['weights'])[1]
    # string += '%.4e	%.4e	' % (value*100, std*100)
    #
    # # sum X
    # value = np.power(10, solution['fit_params']['1H2-16O']['value']) + \
    #         np.power(10, solution['fit_params']['12C-1H4']['value']) + \
    #         np.power(10, solution['fit_params']['14N-1H3']['value']) + \
    #         np.power(10, solution['fit_params']['12C-16O2']['value'])
    # std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']) +
    #                            np.power(10, solution['fit_params']['12C-1H4']['trace']) +
    #                            np.power(10, solution['fit_params']['14N-1H3']['trace']) +
    #                            np.power(10, solution['fit_params']['12C-16O2']['trace']),
    #                            solution['weights'])[1]
    # string += '%.4e	%.4e	' % (value*100, std*100)

    for param in params_others:
        if not param in solution['fit_params']:
            string += 'fix	fix	'
        else:
            if param == 'coupled_mu':
                value = solution['fit_params'][param]['value']/AMU
                std = weighted_avg_and_std(solution['fit_params'][param]['trace'], solution['weights'])[1]/AMU
            else:
                value = solution['fit_params'][param]['value']
                std = weighted_avg_and_std(solution['fit_params'][param]['trace'], solution['weights'])[1]
            string += '%s   %.4e	%.4e	' % (param, value, std)
    print string
