'''
Create pickled sigma array for given molecule from ascii files

Fileformat is exomol v1 or v2 (zero pressure or pressure broadened)

For pressure broadened xsec, temperatures and pressures need to be
consistent. E.g. same temperatures for all pressures and viceversa.

Usage:

python create_xsec.py -s 'source_directory'
                      -o 'output_directory'
                      -b 'binning_resolution'
                      -n 'molecule_name'
                      -v 'exomol_file_version'

'''

import sys, os, optparse, string, pickle, glob
import numpy as np

parser = optparse.OptionParser()
parser.add_option('-s', '--source',
                  dest='source_files',
                  default=None,
)
parser.add_option('-o', '--output',
                  dest='output',
                  default='Input',
)
# binning resolution in wno. If not specified use full resolution
parser.add_option('-b', '--bin',
                  dest='bin',
                  default=None,
)
parser.add_option('-n', '--name',
                  dest='name',
                  default='sigma',
)
# exomol file version. v0: zero pressure, v1: pressure broadened
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

    fname  = os.path.basename(fname)

    if int(filetype) == 1:
        # exomol file format v1 (zero pressure), e.g. 1H2-16O_0-30000_500K_0.010000.sigma
        s = string.split(fname,'_',5)
        temperature = float(s[2][:-1])
        pressure = 0.
        resolution = 0.01 # this should be splitname[5], but need to remove extension

    elif int(filetype) == 2:
        # exomol file format v2 (pressure broadened), e.g. 12C-16O__0-8500_300_0.1_0.01.sig
        s = string.split(fname,'_',5)
        temperature = float(s[3])
        pressure = float(s[4])
        resolution = 0.01 # this should be splitname[5], but need to remove extension

    return temperature, pressure, resolution

def get_specbingrid(wavegrid, specgrid):

    bingrid =[]
    bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
    for i in range(len(wavegrid)-1):
        bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
    bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
    bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins

    return bingrid, bingrid_idx

temperatures = []
pressures = []
wngrid = []

# identify all pressures and temperatures and check for wavenumber grid constitency
for fname in glob.glob(os.path.join(options.source_files, '*.*')):

    temperature, pressure, resolution = read_filename(fname, options.version)
    temperatures.append(temperature)
    pressures.append(pressure)

    print 'Reading %s' % fname

    sigma = np.loadtxt(fname)

    if len(wngrid) == 0:
        wngrid = sigma[:,0]
    else:

        if len(wngrid) <> len(sigma[:,0]):
            print 'Wavenumber grid is not consistent for file %s. Skpping' % os.path.basename(fname)
            continue


print 'Sorting  pressures and temperatures'
pressures = np.sort(np.unique(pressures))
temperatures = np.sort(np.unique(temperatures))
print 'Pressures are %s' % pressures
print 'Temperatures are %s' % temperatures


if bin:
    new_wngrid = np.arange(np.min(wngrid), np.max(wngrid), float(options.bin))
    bin_grid, bin_grid_idx = get_specbingrid(new_wngrid, wngrid)

sigma_array = np.zeros((len(pressures), len(temperatures), len(wngrid)))

for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        for fname in glob.glob(os.path.join(options.source_files, '*.*')):

            xsec_t, xsec_p, xsec_r = read_filename(fname, options.version)
            if xsec_t == temperature_val and xsec_p == pressure_val:
                print fname
                sigma = np.loadtxt(fname)
                if bin:
                    logval = np.log(sigma[:,1])
                    values = np.asarray([np.average(logval[bin_grid_idx == i]) for i in xrange(1,len(bin_grid))])
                    values[np.isnan(values)] = 0
                    values = np.exp(values)
                else:
                    values = sigma[:,1]

                sigma_array[pressure_idx, temperature_idx, :] = sigma[:,1]

sigma_out = {
    'name': options.name,
    'pres': pressures,
    'temp': temperatures,
    'wno': wngrid,
    'sigma': sigma_array,
}

pickle.dump(open(os.path.join(options.output, 'wb'), '%s.db' % options.name), sigma_out)


