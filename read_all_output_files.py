#loading libraries
import pylab, sys, os, optparse, time, pickle
import numpy as np
from ConfigParser import SafeConfigParser
import math

AMU   = 1.660538921e-27 #atomic mass to kg

parser = optparse.OptionParser()
parser.add_option('-i', '--input',
                  dest="grid_folder",
)

options, remainder = parser.parse_args()

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


params_clr = ['1H2-16O_CLR', '12C-1H4_CLR', '14N-1H3_CLR', '12C-16O2_CLR', 'H2_CLR', 'He_CLR', 'N2_CLR']
params_gas = ['1H2-16O', '12C-1H4', '14N-1H3', '12C-16O2', 'H2', 'He', 'N2']
params_others = ['T', 'P0', 'Radius', 'coupled_mu']


print options.grid_folder

for res in os.listdir(options.grid_folder):
    if os.path.isdir(os.path.join(options.grid_folder, res)):
        for error in os.listdir(os.path.join(options.grid_folder, res)):
            if os.path.isdir(os.path.join(options.grid_folder, res, error)):
                filename = os.path.join(options.grid_folder, res, error, 'NEST_out.db')
                if os.path.isfile(filename):

                    out = pickle.load(open(filename, 'rb'))

                    for solution in out:
                        string = '%s	%s	' % (res, error)

                        # for param in params_gas:
                        #     value = solution['fit_params'][param]['value']
                        #     std = weighted_avg_and_std(solution['fit_params'][param]['trace'], solution['weights'])
                        #     string += '%.4e	%.4e	' % (value, std)

                        for param in params_gas:
                            value = np.power(10, solution['fit_params'][param]['value'])
                            std = weighted_avg_and_std(np.power(10,solution['fit_params'][param]['trace']), solution['weights'])[1]
                            string += '%.4e	%.4e	' % (value*100, std*100)

                        # H20/CH4
                        value = np.power(10, solution['fit_params']['1H2-16O']['value']-solution['fit_params']['12C-1H4']['value'])
                        std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']-solution['fit_params']['12C-1H4']['trace']),
                                                   solution['weights'])[1]
                        string += '%.4e	%.4e	' % (value, std)

                        # H2O/NH3
                        value = np.power(10, solution['fit_params']['1H2-16O']['value']-solution['fit_params']['14N-1H3']['value'])
                        std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']-solution['fit_params']['14N-1H3']['trace']),
                                                   solution['weights'])[1]
                        string += '%.4e	%.4e	' % (value, std)

                        # CH4/NH3
                        value = np.power(10, solution['fit_params']['12C-1H4']['value']-solution['fit_params']['14N-1H3']['value'])
                        std = weighted_avg_and_std(np.power(10, solution['fit_params']['12C-1H4']['trace']-solution['fit_params']['14N-1H3']['trace']),
                                                   solution['weights'])[1]
                        string += '%.4e	%.4e	' % (value, std)

                        #He+H2+N2
                        value = np.power(10, solution['fit_params']['He']['value'])+ np.power(10, solution['fit_params']['H2']['value']) + \
                                np.power(10, solution['fit_params']['N2']['value'])
                        std = weighted_avg_and_std(np.power(10, solution['fit_params']['He']['trace']) +
                                                   np.power(10, solution['fit_params']['H2']['trace']) +
                                                   np.power(10, solution['fit_params']['N2']['trace']),
                                                   solution['weights'])[1]
                        string += '%.4e	%.4e	' % (value*100, std*100)

                        # sum X
                        value = np.power(10, solution['fit_params']['1H2-16O']['value']) + \
                                np.power(10, solution['fit_params']['12C-1H4']['value']) + \
                                np.power(10, solution['fit_params']['14N-1H3']['value']) + \
                                np.power(10, solution['fit_params']['12C-16O2']['value'])
                        std = weighted_avg_and_std(np.power(10, solution['fit_params']['1H2-16O']['trace']) +
                                                   np.power(10, solution['fit_params']['12C-1H4']['trace']) +
                                                   np.power(10, solution['fit_params']['14N-1H3']['trace']) +
                                                   np.power(10, solution['fit_params']['12C-16O2']['trace']),
                                                   solution['weights'])[1]
                        string += '%.4e	%.4e	' % (value*100, std*100)

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
                                string += '%.4e	%.4e	' % (value, std)
                        print string
