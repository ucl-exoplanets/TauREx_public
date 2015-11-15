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
                  default='xsec',
                  )
parser.add_option('-e', '--extension',
                  dest='extension',
                  default='sigma',
                  )
parser.add_option('-l', '--wnlim',
                  dest='wnlim',
                  default=None,
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

    bingrid = np.asarray([wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0 for i in range(len(wavegrid)-1)])
    bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
    return bingrid, bingrid_idx

temperatures = []
pressures = []
wngrid = []

# identify all pressures and temperatures and check for wavenumber grid constitency
for fname in glob.glob(os.path.join(options.source_files, '*.%s' % options.extension)):

    temperature, pressure, resolution = read_filename(fname, options.version)
    temperatures.append(temperature)
    pressures.append(pressure)

    print 'Reading %s' % fname

    sigma = np.loadtxt(fname)

    if len(wngrid) == 0:
        wngrid = sigma[:,0]
        if options.bin == None:
            out_wngrid = wngrid
            if options.wnlim:
                out_wngrid = out_wngrid[np.logical_and((out_wngrid > options.wnmin[0]), out_wngrid<options.wnmax[0])]
    else:
        if len(wngrid) <> len(sigma[:,0]):
            print 'Wavenumber grid is not consistent for file %s. Skpping' % os.path.basename(fname)
            continue

if options.wnlim:
    wnlim = options.wnlim.split(',')
    wnmin = float(wnlim[0])
    wnmax = float(wnlim[1])
    wngrid = wngrid[np.logical_and(wngrid>wnmin, wngrid<wnmax)]
else:
    wnmin = np.min(wngrid)
    wnmax = np.max(wngrid)

print 'Sorting  pressures and temperatures'
pressures = np.sort(np.unique(pressures))
temperatures = np.sort(np.unique(temperatures))
print 'Pressures are %s' % pressures
print 'Temperatures are %s' % temperatures

if options.bin:
    new_wngrid = np.arange(wnmin, wnmax, float(options.bin))
    bin_grid, bin_grid_idx = get_specbingrid(new_wngrid, wngrid)
    out_wngrid = bin_grid
    if options.wnlim:
        out_wngrid = out_wngrid[np.logical_and(out_wngrid>wnmin, out_wngrid<wnmax)]


sigma_array = np.zeros((len(pressures), len(temperatures), len(bin_grid)))

for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        for fname in glob.glob(os.path.join(options.source_files, '*.%s' % options.extension)):
            xsec_t, xsec_p, xsec_r = read_filename(fname, options.version)
            if xsec_t == temperature_val and xsec_p == pressure_val:
                print 'Binning in log space (geometric average): %s' % fname
                sigma = np.loadtxt(fname)
                if options.bin:
                    logval = np.log(sigma[:,1])
                    values = np.asarray([np.average(logval[bin_grid_idx == i]) for i in xrange(0,len(bin_grid))])
                    values[np.isnan(values)] = 0
                    values = np.exp(values)
                else:
                    values = sigma[:,1]
                if options.wnlim:
                    values = values[np.logical_and(bin_grid>wnmin, bin_grid<wnmax)]
                sigma_array[pressure_idx, temperature_idx, :] = values

sigma_out = {
    'name': options.name,
    'p': pressures,
    't': temperatures,
    'wno': out_wngrid,
    'xsecarr': sigma_array,
}

pickle.dump(sigma_out, open(os.path.join(options.output, '%s.db' % options.name), 'wb'))