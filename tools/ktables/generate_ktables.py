'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    Generate k tables

    Developers: Marco Rocchetto (University College London)

'''

import numpy as np
import cPickle as pickle
import os
from ConfigParser import SafeConfigParser
import argparse
import time
import glob
import sys

# support functions
def get_specgrid( R=5000, lambda_min=0.1, lambda_max=20.0):
    #generating wavelength grid with uniform binning in log(lambda)
    #lambda min and max are in microns, R is the spectral resolution
    #R = lambda/delta_lamdba
    specgrid = []
    delta_lambda =[]
    specgrid.append(lambda_min)
    run = True
    i=1
    while run:
        dlam= specgrid[i-1]/R
        specgrid.append(specgrid[i-1]+dlam)
        delta_lambda.append(dlam)

        if specgrid[i] >= lambda_max:
            run=False
        i+=1
    return np.asarray(specgrid),np.asarray(delta_lambda)


#loading parameter file parser
parser = argparse.ArgumentParser()


parser.add_argument('--input_file',
                      dest='input_file'
                    )

parser.add_argument('--input_files',
                      dest='input_files',
                      default=False,
                    )

parser.add_argument('--output_file',
                      dest='output_filename',
                      default=False,
                    )

parser.add_argument('--input_ext',
                      dest='input_ext',
                      default=False,
                    )
parser.add_argument('--output_folder',  # can be pickle or ascii
                      dest='output_folder',
                      default=False
                    )

parser.add_argument('--wl',
                      dest='wlrange',
                      default=False,
                    )
parser.add_argument('--wn',
                      dest='wnrange',
                      default=False,
                    )
parser.add_argument('--resolution',         # if wnrange is set, set the wavenumber bin (will be constant)
                      dest='resolution',    # if wlrange is set, set the resolution in wavelength space (will be constant)
                      default=False,
                    )

parser.add_argument('--ngauss',
                      dest='ngauss',
                      default=20,
                    )

parser.add_argument('--p_idx',
                      dest='p_idx',
                      default=False,
                    )
parser.add_argument('--t_idx',
                      dest='t_idx',
                      default=False,
                    )
parser.add_argument('--cores',
                      dest='cores',
                      default=1,
                    )

options = parser.parse_args()

if not options.input_file and not options.input_files:
    print 'You need to specify an input file or folder'
    exit()
if options.input_files and not options.output_folder:
    print 'You need to specify an output folder for the ktables'
    exit()
if not options.wlrange and not options.wnrange:
    print 'You need to specify an input xsec'
    exit()
if  options.wlrange and options.wnrange:
    print 'You need to set either --wl (wavelength range) or --wn (wavenumber range)'
    exit()
if not options.resolution:
    print 'You need to specify a resolution'
    exit()


start = startall = time.time()

if options.wlrange:
    lambda_min, lambda_max = options.wlrange.split(',')
    wl_grid = get_specgrid(int(options.resolution), float(lambda_min), float(lambda_max))[0] # build bin grid in wavelength
    wn_grid = np.sort(10000./wl_grid) # convert to wavenumber and sort
    bincentres = np.sort(10000./(0.5*(wl_grid[1:] + wl_grid[:-1]))) # get the bin centres

elif options.wnrange:
    wn_min, wn_max = options.wnrange.split(',')
    wn_grid = np.arange(float(wn_min), float(wn_max), float(options.resolution))
    bincentres = 0.5*(wn_grid[1:] + wn_grid[:-1]) # get the bin centres
ngauss = int(options.ngauss)
gauss = np.polynomial.legendre.leggauss(ngauss) # get the legendre-gauss sample points and weights
samples = gauss[0]
weights = gauss[1]
end = time.time()

files = []

if options.input_file:
    files.append(options.input_file)
else:
    files = glob.glob(os.path.join(options.input_files))

start = time.time()

for filenm in files:

    try:
        xsec = pickle.load(open(filenm))
    except:
        xsec = np.loadtxt(filenm)

    print 'Process %s...' % filenm,
    try:
        bingrid_idx = np.digitize(xsec[:,0] ,wn_grid) # get the indexes for bins

        # build the k table
        ktable = np.zeros((len(wn_grid)-1, ngauss))

        for i in range(1, len(wn_grid)):
            x = xsec[:,0][bingrid_idx == i]
            abs = xsec[:,1][bingrid_idx == i]
            sort_bin = np.sort(abs)
            if len(np.unique(sort_bin)) > 1:
                norm_x = ((x-np.min(x))/(np.max(x) - np.min(x)))*2 - 1
                kcoeff = np.interp(gauss[0], norm_x, sort_bin)
                ktable[i-1,:]  = kcoeff

        kdist_out = {}
        kdist_out['wno'] = bincentres
        kdist_out['resolution'] = options.resolution
        kdist_out['weights'] = weights
        kdist_out['ngauss'] = ngauss
        kdist_out['kcoeff'] = ktable

        if options.t_idx or options.p_idx:
            split = os.path.splitext(os.path.basename(filenm))[0].split('_')
        if options.p_idx:
            try:
                kdist_out['p'] = float(split[int(options.p_idx)])
            except:
                try:
                    kdist_out['p'] = float(split[int(options.p_idx)][1:])
                except:
                    pass
        if options.t_idx:
            try:
                print split[int(options.t_idx)]
                kdist_out['t'] = float(split[int(options.t_idx)])
            except:
                try:
                    print split[int(options.t_idx)][1:]
                    kdist_out['t'] = float(split[int(options.t_idx)][1:])
                except:
                    pass

        end = time.time()
        print 'Done in %.3f s' % (end-start)

        if options.output_filename and options.input_file and not options.input_files:
            out_filename = options.output_filename
        else:
            out_filename = '%s.ktable' % os.path.basename(os.path.splitext(filenm)[0])

        if options.output_folder:
            out_filename = os.path.join(os.path.abspath(options.output_folder), out_filename)
        elif not options.input_files:
            out_filename = os.path.join(os.path.dirname(os.path.abspath(options.input_file)), out_filename)

        start = time.time()
        print 'Dump ktable to %s...' % out_filename,
        pickle.dump(kdist_out, open(out_filename, 'wb'), protocol=2)
        end = time.time()
        print 'Done in %.3f s' % (end-start)
    except:
        print "Unexpected error:", sys.exc_info()
        pass

endall = time.time()

print 'Stats: %i ktables computed in %.3f s' % (len(files), (endall-startall))