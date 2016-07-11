'''

Bin cross sections

'''

import sys, os, optparse, string, glob
import numpy as np
from time import gmtime, strftime
import multiprocessing
from multiprocessing import Pool
from functools import partial


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
parser.add_option('-c', '--cores',
                  dest='cores',
                  default=1,
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


def fast_nearest_interp(xi, x, y):
    """Assumes that x is monotonically increasing!!."""
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

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


def round_base(x, base=.05):
  return base * round(float(x)/base)



files = glob.glob(os.path.join(options.source_files, '*.%s' % options.extension))

print 'Create wavenumber binning grid'

fname = files[0]
temperature, pressure, res_file = read_filename(fname, options.version)
wngrid = np.loadtxt(fname)[:,0]
if options.linear_binning:
    # create the output wavenumber grid
    bin_wngrid = np.arange(round_base(np.min(wngrid), float(options.linear_binning)),
                           round_base(np.max(wngrid), float(options.linear_binning)),
                           float(options.linear_binning))
    bingrid_idx = np.digitize(wngrid, bin_wngrid) #getting the specgrid indexes for bins
    wngrid = bin_wngrid

print 'Done'

def worker(idx):

    files = glob.glob(os.path.join(options.source_files, '*.%s' % options.extension))

    fname = files[idx]

    xsec_t, xsec_p, xsec_r = read_filename(fname, options.version)

    filename = os.path.join(options.output, 'binned_P%.2e_T%i_%s_%.2f.xsec' % (xsec_p,
                                                                               xsec_t,
                                                                               options.molecule_name,
                                                                               float(options.linear_binning)))
    if os.path.isfile(filename):
        print 'Skip file %i/%i: %s' %  (idx+1, len(files), fname)
    else:
        print 'Bin file %i/%i: %s' %  (idx+1, len(files), fname)

        loadf = np.loadtxt(fname)
        sigma_in = loadf[:,1]

        if options.binning_method == 'geometric_average':
            # geometric average (i.e. log-average)
            logval = np.log(sigma_in)
            values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid)+1)])
            values = np.exp(values)
            values[np.isnan(values)] = 0
        elif options.binning_method == 'algebraic_average':
            # algebraic average
            values = np.asarray([np.average(sigma_in[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid)+1)])
            values[np.isnan(values)] = 0
        elif options.binning_method == 'random_sample':
            values = np.asarray([np.random.choice(sigma_in[bingrid_idx == i], 1)[0] for i in xrange(1,len(bin_wngrid)+1)])
        elif options.binning_method == 'first_sample':
            values = np.asarray([sigma_in[bingrid_idx == i][0] for i in xrange(1,len(bin_wngrid)+1)])
        elif options.binning_method == 'interp_sample':
            values = np.interp(bin_wngrid, loadf[:,0], sigma_in)
        elif options.binning_method == 'interp_nearest_sample':
            values = fast_nearest_interp(bin_wngrid, loadf[:,0], sigma_in)

        out = np.zeros((len(values), 2))
        out[:,0] = wngrid
        out[:,1] = values
        np.savetxt(filename, out)



pool = Pool(processes=int(options.cores))      #setting number of cores on which to run
pool_result = pool.map(worker,[i for i in range(len(files))]) #runnning the stuff
pool.close()                            #closing the pool 
pool.join()                             #closing the pool. Needs to be done otherwise createob.reset has no effect.

# 
# nrep = len(files)/int(options.cores)
# 
# procs = []
# count = 0
# for i in range(nrep):
#     for n in range(int(options.cores)):
#         p = multiprocessing.Process(target=worker, args=(count, ))
#         p.start()
#         procs.append(p)
#         count += 1
#     for p in procs:
#         p.join()
#     procs = []
# 
# for i in range(len(files) % int(options.cores)):
#     for n in range(int(options.cores)):
#         p = multiprocessing.Process(target=worker, args=(count, ))
#         p.start()
#         count += 1
#     for p in procs:
#         p.join()
#     procs = []
