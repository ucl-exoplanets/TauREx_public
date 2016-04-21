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
                      -o 'output_filename'
                      -n 'molecule_name'
                      -v 'exomol_file_version' (1: zero pressure, 2: pressure broadened)
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
parser.add_option('-o', '--output',
                  dest='output',
                  default=None,
)
parser.add_option('-n', '--molecule_name',
                  dest='molecule_name',
                  default=None,
                  )
parser.add_option('-e', '--extension',
                  dest='extension',
                  default='sigma',
                  )
parser.add_option('-b', '--binning_resolution',
                  dest='linear_binning',
                  default=None,
)
parser.add_option('-m', '--binning_method',
                  dest='binning_method',
                  default='geometric_average',
)
# exomol file version. ZP: zero pressure, PB: pressure broadened, CH4: Sergey's ch4
parser.add_option('-v', '--version',
                  dest='version',
                  default='ZP',
)
options, remainder = parser.parse_args()

if not options.source_files:
    print 'Select source files with option -s'
    exit()
if not options.linear_binning:
    print 'Select a binning resolution                                        '
    exit()

# compatibility with exomol filename formats
def read_filename(fname, filetype=1):

    fname  = os.path.splitext(os.path.basename(fname))[0]

    if filetype.upper() == 'ZP':
        # exomol file format Zero Pressure, e.g. 1H2-16O_0-30000_500K_0.010000.sigma
        s = string.split(fname,'_',5)
        temperature = float(s[2][:-1])
        pressure = 0.
        resolution = float(s[3])

    elif filetype.upper() == 'PB':
        # exomol file format Pressure Broadened, e.g. 12C-16O__0-8500_300_0.1_0.01.sig
        strfname = fname.replace('__', '_')
        s = strfname.split('_')
        temperature = float(s[2])
        pressure = float(s[3])
        resolution = float(s[4])

    elif filetype.upper() == 'EX':
        s = fname.split('_')
        pressure = float(s[1][1:])
        temperature = float(s[2][1:])
        resolution = 0.01 # assume this resolution...

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
print 'Wavenumber range is %f - %f' % (np.min(wngrid), np.max(wngrid))

comments.append('The resolution of the original Exomol cross sections is %f' % resolution)

def round_base(x, base=.05):
  return base * round(float(x)/base)

if options.linear_binning:
    # create the output wavenumber grid
    bin_wngrid = np.arange(round_base(np.min(wngrid), float(options.linear_binning)),
                           round_base(np.max(wngrid), float(options.linear_binning)),
                           float(options.linear_binning))
    bingrid_idx = np.digitize(wngrid, bin_wngrid) #getting the specgrid indexes for bins

    wngrid = bin_wngrid

    comments.append('Linear binning, at %f wavenumber resolution.' % float(options.linear_binning))

    if options.binning_method == 'geometric_average':
        comments.append('Linear binning: use geometric average.' )
    elif options.binning_method == 'algebraic_average':
        comments.append('Linear binning: use algebraic average.' )

for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        for fname in glob.glob(os.path.join(options.source_files, '*.%s' % options.extension)):

            xsec_t, xsec_p, xsec_r = read_filename(fname, options.version)
            if xsec_t == temperature_val and xsec_p == pressure_val:

                filename = os.path.join(options.output, 'binned_P%.2e_T%i_%s_%.2f.xsec' % (pressure_val,
                                                                                           temperature_val,
                                                                                           options.molecule_name,
                                                                                           float(options.linear_binning)))
                if os.path.isfile(filename):
                    print 'Skip file %i/%i: %s' % (pressure_idx*temperature_idx,
                                                  len(pressures)*len(temperatures),
                                                  fname)
                else:
                    print 'Bin file %i/%i: %s' % (pressure_idx*temperature_idx,
                                                  len(pressures)*len(temperatures),
                                                  fname)
                    sigma_in = np.loadtxt(fname)[:,1]

                    if options.binning_method == 'geometric_average':
                        # geometric average (i.e. log-average)
                        logval = np.log(sigma_in)
                        values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid))+1])
                        values = np.exp(values)
                        values[np.isnan(values)] = 0
                    elif options.binning_method == 'algebraic_average':
                        # algebraic average
                        values = np.asarray([np.average(sigma[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid)+1)])
                        values[np.isnan(values)] = 0

                    out = np.zeros((len(values), 2))
                    out[:,0] = wngrid
                    out[:,1] = values
                    np.savetxt(filename, out)



