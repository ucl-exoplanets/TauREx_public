
# python create_hitran_xsec.py -o '/Volumes/DATA/xsec/pressure_broadened_xsec/hitran/C2H2/C2H2_hitran_PB.db' -n C2H2 -w 605,6889 -r 0.01


from hapi import *

import sys, os, optparse, string, pickle, glob
import numpy as np
from time import gmtime, strftime

parser = optparse.OptionParser()

parser.add_option('-o', '--output_filename',
                  dest='output_filename',
                  default='sigma_array.db',
)

parser.add_option('-n', '--molecule_name',
                  dest='molecule_name',
                  default=None,
                  )
parser.add_option('-w', '--wavenumber_range',
                  dest='wavenumber_range',
                  default=None,
                  )
parser.add_option('-r', '--resolution',
                  dest='resolution',
                  default=0.01,
)

options, remainder = parser.parse_args()

wnmin, wnmax = options.wavenumber_range.split(',')
wnmin = float(wnmin)
wnmax = float(wnmax)

p_bar = np.asarray([1.00000000e-03,   3.00000000e-03,   5.00000000e-03,
         1.00000000e-02,   2.00000000e-02,   4.00000000e-02,
         8.00000000e-02,   1.00000000e-01,   3.00000000e-01,
         6.00000000e-01,   9.00000000e-01,   1.00000000e+00,
         1.20000000e+00,   1.50000000e+00,   2.00000000e+00,
         2.50000000e+00,   3.00000000e+00,   4.00000000e+00,
         5.00000000e+00,   6.00000000e+00,   8.00000000e+00,
         1.00000000e+01])

temperatures =  np.asarray([  300.,   400.,   500.,   600.,   700.,   800.,   900.,  1000.,
        1100.,  1200.,  1300.,  1400.,  1500.])

pressures = p_bar * 0.986923

fetch(options.molecule_name, 1, 1, wnmin, wnmax)


# get wavenumber grid
wngrid, coeff = absorptionCoefficient_Voigt(SourceTables=options.molecule_name, Environment={'T':temperatures[0],'p':pressures[0]}, OmegaWing=10, OmegaStep=float(options.resolution))

sigma_out = {
    'name': options.molecule_name,
    'p': p_bar,
    't': temperatures,
    'wno': wngrid,
}

sigma_array = np.zeros((len(pressures), len(temperatures), len(wngrid)))

for pres_idx, pres_val in enumerate(pressures):
    for temp_idx, temp_val in enumerate(temperatures):

        print 'Computing cross section for T = %i, P = %.4f' % (temp_val, pres_val)

        nu, coeff = absorptionCoefficient_Voigt(SourceTables=options.molecule_name, Environment={'T':temp_val,'p':pres_val})
        sigma_array[pres_idx, temp_idx, :] = coeff

sigma_out['xsecarr'] = sigma_array

comments = []
comments.append('Cross section array created from HITRAN cross sections, using create_hitran_xsec.py and HAPI, on GMT %s' %
                (strftime("%Y-%m-%d %H:%M:%S", gmtime())))
comments.append('The resolution of the HITRAN cross sections is 0.01 wn')

sigma_out['comments'] = comments

pickle.dump(sigma_out, open(options.output_filename, 'wb'))

