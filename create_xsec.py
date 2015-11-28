'''

Create pickled high resolution sigma array for given molecule from ascii files

Fileformat is exomol v1 or v2 (zero pressure or pressure broadened).

NB:

- Use the highest resolution available (<=0.01)

- The pressure and temperature grids are the ones provided by the input files.

- For pressure broadened xsec, temperatures and pressures need to be
  consistent. E.g. same temperatures for all pressures and viceversa.

- If you want to reinterpolate to a different pressure/temperature grid, firstly
  use this script to convert exomol xsec to a pickled xsec, then use convert_xsec.py.

- Similarly, if you want to bin these high res xsec to a different wavenumber grid
  (linear or nonlinear), firstly use this script, then use convert_xsec.py.


Usage:

python create_xsec.py -s 'source_directory'
                      -o 'output_directory'
                      -n 'molecule_name'
                      -v 'exomol_file_version'
                      -e 'file_extension'
'''

import sys, os, optparse, string, pickle, glob
import numpy as np
from time import gmtime, strftime

parser = optparse.OptionParser()
parser.add_option('-s', '--source',
                  dest='source_files',
                  default=None,
)
parser.add_option('-o', '--output_filename',
                  dest='output_filename',
                  default='sigma_array.db',
)
# binning resolution in wno. If not specified use full resolution
# parser.add_option('-b', '--bin',
#                   dest='bin',
#                   default=None,
# )
parser.add_option('-n', '--molecule_name',
                  dest='molecule_name',
                  default=None,
                  )
parser.add_option('-e', '--extension',
                  dest='extension',
                  default='sigma',
                  )

# exomol file version. 1: zero pressure, 2: pressure broadened
parser.add_option('-v', '--version',
                  dest='version',
                  default=1,
)
options, remainder = parser.parse_args()

if not options.source_files:
    print 'Select source files with option -s'
    exit()

# compatibility with exomol filename formats
def read_filename(fname, filetype=1):

    fname  = os.path.splitext(os.path.basename(fname))[0]

    if int(filetype) == 1:
        # exomol file format v1 (zero pressure), e.g. 1H2-16O_0-30000_500K_0.010000.sigma
        s = string.split(fname,'_',5)
        temperature = float(s[2][:-1])
        pressure = 0.
        resolution = s[3]

    elif int(filetype) == 2:
        # exomol file format v2 (pressure broadened), e.g. 12C-16O__0-8500_300_0.1_0.01.sig
        s = string.split(fname,'_',5)
        temperature = float(s[3])
        pressure = float(s[4])
        resolution = s[5]

    return temperature, pressure, resolution

resolution = None
comments = []
temperatures = []
pressures = []
wngrid = []


# identify all pressures and temperatures and check for wavenumber grid constitency
for fname in glob.glob(os.path.join(options.source_files, '*.%s' % options.extension)):

    temperature, pressure, res_file = read_filename(fname, options.version)
    temperatures.append(temperature)
    pressures.append(pressure)

    if not resolution:
        resolution = res_file
    else:
        if res_file != resolution:
            print 'Resolution is not consistent for file %s. Skpping' % os.path.basename(fname)
            continue

    print 'Reading %s' % fname

    if len(wngrid) == 0:
        sigma = np.loadtxt(fname)
        wngrid = sigma[:,0]


print 'Sorting  pressures and temperatures'
pressures = np.sort(np.unique(pressures))
temperatures = np.sort(np.unique(temperatures))
print 'Pressures are %s' % pressures
print 'Temperatures are %s' % temperatures

sigma_array = np.zeros((len(pressures), len(temperatures), len(wngrid)))

for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        for fname in glob.glob(os.path.join(options.source_files, '*.%s' % options.extension)):
            xsec_t, xsec_p, xsec_r = read_filename(fname, options.version)
            if xsec_t == temperature_val and xsec_p == pressure_val:
                print 'Reading %s' % fname
                values = np.loadtxt(fname)[:,1]
                sigma_array[pressure_idx, temperature_idx, :] = values

time = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
comments.append('Cross section array created from Exomol cross sections v%i, using create_xsec.py, on GMT %s' %
                (options.version, time))
comments.append('The resolution of the Exomol cross sections is %f' % resolution)

sigma_out = {
    'name': options.name,
    'p': pressures,
    't': temperatures,
    'wno': out_wngrid,
    'xsecarr': sigma_array,
    'comments': comments
}

pickle.dump(sigma_out, open(options.output_filename), 'wb')