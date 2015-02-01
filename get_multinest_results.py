################################################
# get_multinest_results.py
#
# This scripts creates a text file with output results from
# a multinest run.
#
# Input: Output folder
#
# Output: The scripts writes into a output_multinest.txt, placed inside the output folder
#
# Modification History:
#   - v1.0 : first definition, Marco Rocchetto 31 Jan 2014
#
################################################



import numpy as np
import pymultinest
import optparse
import os

parser = optparse.OptionParser()
parser.add_option('-o', '--output',
                  dest="output",
                  default="./Output",
)
options, remainder = parser.parse_args()

# parameters.dat contains the list of fitted parameters
fit_params = np.loadtxt(os.path.join(options.output, 'parameters.dat') ,dtype='string')

# get multinest data & stats
analyzer = pymultinest.Analyzer(n_params=len(fit_params),
                                outputfiles_basename=os.path.join(options.output, 'multinest', '1-'))
stats = analyzer.get_stats()

# create list of output results for the different modes detected
out_list = []
out_list_std = []
if len(stats['modes']) > 1:
    # more than one mode found
    for n in range(len(stats['modes'])):
        mode_list = []
        mode_list_std = []
        for i in range(len(fit_params)):
            mode_list.append(stats['modes'][n]['maximum a posterior'][i])
            mode_list_std.append(stats['modes'][n]['sigma'][i])
        out_list.append(mode_list)
        out_list_std.append(mode_list_std)

# inverse transformation of clr(X) [mixing ratios are in the cenetered-log-ratio space]
params_nogas = ['Temperature', 'Radius', 'P0', 'mu'] # all parameter names excluding gases

out_mixing = []
out_mixing_std = []
for idx_mode, mode in enumerate(out_list):
    out_mixing_mode = []
    out_mixing_std_mode = []
    for idx_param, param in enumerate(mode):
        # build list of mixing ratios
        if not fit_params[idx_param] in params_nogas: #it's a gas ratio
            out_mixing_mode.append(param)
            out_mixing_std_mode.append(out_list_std[idx_mode][idx_param])

    # transform back to standard simplex space
    out_mixing_mode = np.exp(np.asarray(out_mixing_mode))
    out_mixing_mode = out_mixing_mode/np.sum(out_mixing_mode)
    out_mixing_std_mode = out_mixing_mode*np.asarray(out_mixing_std_mode)
    out_mixing.append(out_mixing_mode)
    out_mixing_std.append(out_mixing_std_mode)

# write out to file
with open(os.path.join(options.output, 'output_multinest.txt'), 'w') as file:
    for idx_mode, mode in enumerate(out_list):
        file.write('Mode %i: \n' % idx_mode)

        ngas = 0
        for idx_param, param in enumerate(mode):

            if fit_params[idx_param] in params_nogas: # it's not a gas ratio
                file.write('%s: %.12f +/- %.12f  \n' % (fit_params[idx_param], param, out_list_std[idx_mode][idx_param]))
            else:  #it's a gas ratio
                file.write('%s: %.12f +/- %.12f  \n' % (fit_params[idx_param],
                                                      out_mixing[idx_mode][ngas],
                                                      out_mixing_std[idx_mode][ngas]))
                ngas += 1
        file.write('\n')




























