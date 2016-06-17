'''
Reinterpolate a pickled sigma array created with create_xsec.py or otherwise
to a user defined temperaturem, pressure and wavenumber grid.



It is also possible to bin the cross section (geometric average method), with
constant user defined bin sizes.

NB:

- Always use the highest resolution input cross section (at least 0.01 wn resolution!)
  to obtain reliable results...

- The bins of the input must be constant (e.g. 0.01 wavenumbers), and specified either
  in the sigma_array (sigma_array['resolution'], or with the optional parameter '-r'

- if the requested temperature or pressure grids extends over/below the temperature or
  pressure ranges available in the input pickled xsec, the highest/lowest temperature or
  pressure available in the input xsec will be used.



Usage:

python convert_xsec.py -s 'sourcce_sigma'
                       -o 'output_filename'
                       -n 'molecule_name'
                       -p 'pressure_list' (comma separated)
                       -t 'temperature_list' (comma separated)
                       -w 'wavenumber_range' (comma separated) [optional]
                       -b 'binning_resolution' [optional]
                       -m 'binning_method' (available are: geometric_average/algebric_average)
                       -r 'input_resolution' [optional]

'''

import sys, os, argparse, string, pickle, glob
import numpy as np
#from scipy.inteprolate import interp2d
from time import gmtime, strftime

import matplotlib.pylab as plt

import multiprocessing

### support functions
def binspectrum(spectrum_in, resolution):
    wavegrid, dlamb_grid = get_specgrid(R=resolution,lambda_min=np.min(spectrum_in[:,0]),lambda_max=np.max(spectrum_in[:,0]))
    spec_bin_grid, spec_bin_grid_idx = get_specbingrid(wavegrid, spectrum_in[:,0])
    spectrum_binned = [spectrum_in[:,1][spec_bin_grid_idx == i].mean() for i in xrange(1,len(spec_bin_grid))]
    return transpose(vstack((wavegrid, spectrum_binned)))

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

def get_specbingrid(wavegrid, specgrid, binwidths=None):
    #function calculating the bin boundaries for the data
    #this is used to bin the internal spectrum to the data in fitting module

    if not isinstance(binwidths, (np.ndarray, np.generic)):
        bingrid =[]
        bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
        for i in range(len(wavegrid)-1):
            bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
        bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
        bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
    else:

        # this bingrid is actually useless, as it doesn't allow for gaps in the data
        bingrid = []
        for i in range(len(wavegrid)):
            bingrid.append(wavegrid[i]-binwidths[i]/2.)
            bingrid.append(wavegrid[i]+binwidths[i]/2.)

        # build bin grid index array (an index for each model datapoint. If the point is outside the input
        # spectrum bins, the idx is -99999 (invalid integer...)
        bingrid_idx = np.empty(len(specgrid))
        bingrid_idx[:] = np.NaN

        for i in range(len(specgrid)):
            for j in range(len(wavegrid)):
                if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                    bingrid_idx[i] = j+1
                    break

    return bingrid, bingrid_idx




parser = argparse.ArgumentParser()
parser.add_argument('-n', '--molecule_name',
                  dest='molecule_name',
                  default=None,
)
parser.add_argument('-s', '--source_sigma',
                  dest='source_sigma',
                  default=None,
)
parser.add_argument('-o', '--output_filename',
                  dest='output_filename',
                  default='sigma_array.db',
)
parser.add_argument('-p', '--pressure_list',
                  dest='pressure_list',
                  default=None,
                  )
parser.add_argument('-t', '--temperature_list',
                  dest='temperature_list',
                  default=None,
                  )
parser.add_argument('-w', '--wavenumber_range',
                  dest='wavenumber_range',
                  default=None,
                  )
parser.add_argument('-l', '--lambda_range',
                  dest='lambda_range',
                  default=None,
                  )

parser.add_argument('--binning_method',
                  dest='binning_method',
                  default='const_wn', # const_wn (constant wavenumber spacing) or const_res (constant resolution in wavelength)
)
parser.add_argument('--binning_resolution',
                  dest='binning_resolution',
                  default=None,
)
parser.add_argument('--sampling_method',
                  dest='sampling_method',
                  default='random_sample',
)
parser.add_argument('--input_resolution',
                  dest='input_resolution',
                  default=None,
)
options = parser.parse_args()

if not options.source_sigma or \
   not options.output_filename:
       print 'Wrong input. Retry...'
       exit()


sigma_in = pickle.load(open(options.source_sigma))

if not options.pressure_list:
   pressures = sigma_in['p']
else:
   l = options.pressure_list.split(',')
   pressures = np.asarray([float(m) for m in l])

if not options.temperature_list:
   temperatures = sigma_in['t']
else:
    l = options.temperature_list.split(',')
    temperatures = np.asarray([float(m) for m in l])

pressures = np.sort(pressures)
temperatures = np.sort(temperatures)

# define the wavenumber range of the output cross sections
# if these limits extend above/below the input xsec, set to 0
if not options.wavenumber_range:
    wmin = np.min(sigma_in['wno'])
    wmax = np.max(sigma_in['wno'])
else:
    wnmin, wnmax = options.wavenumber_range.split(',')
    wnmin = float(wnmin)
    wnmax = float(wnmax)

if options.lambda_range:
    lambdamin, lambdamax = options.lambda_range.split(',')
    lambdamin = float(lambdamin)
    lambdamax = float(lambdamax)


if not options.input_resolution:
    print 'You must specify the resolution in wavenumber of the input sigma array'
else:
    resolution = float(options.input_resolution)

if 'comments' in sigma_in:
    comments = sigma_in['comments']
else:
    comments = []

comments.append('Cross section array created using convert_xsec.py on GMT %s' % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
comments.append('Original wavenumber grid limits were: %f - %f' % (np.min(sigma_in['wno']), np.max(sigma_in['wno'])))
comments.append('Original pressures were: %s' % sigma_in['p'])
comments.append('Original temperatures were: %s' % sigma_in['t'])

print 'Reinterpolate to the requested wavenumber grid from wnmin to wnamx at the input resolution'

# reinterpolate to the requested wavenumber grid from wnmin to wnamx at the input resolution.
# this is done to expand, or restrict, the wavenumber range. Zeros are assumed if we extend the wavenumber range.
full_wngrid = np.arange(wnmin, wnmax, resolution)
sigma_array = np.zeros((len(sigma_in['p']), len(sigma_in['t']), len(full_wngrid)))
for pressure_idx, pressure_val in enumerate(sigma_in['p']):
    for temperature_idx, temperature_val in enumerate(sigma_in['t']):
        sigma_array[pressure_idx, temperature_idx] =  np.interp(full_wngrid, sigma_in['wno'], sigma_in['xsecarr'][pressure_idx, temperature_idx], left=0, right=0)


# bin if requested

if options.binning_method == 'const_res':

    print 'Constant resolution in lambda R=%.2f binning' % float(options.binning_resolution)
    comments.append('Constant resolution in lambda R=%.2f binning' % float(options.binning_resolution))

    if options.sampling_method == 'geometric_average':
        comments.append('Linear binning: use geometric average.' )
    elif options.sampling_method == 'algebraic_average':
        comments.append('Linear binning: use algebraic average.' )
    elif options.sampling_method == 'random_sample':
        comments.append('Linear binning: use random sample.' )
    elif options.sampling_method == 'first_sample':
        comments.append('Linear binning: use first sample.' )

    if wnmin == 0:
        wnmin = 1

    wavegrid, dlamb_grid = get_specgrid(R=float(options.binning_resolution),lambda_min=lambdamin,lambda_max=lambdamax)

    if options.sampling_method == 'geometric_average' or options.sampling_method == 'algebraic_average':
        spec_bin_grid, bingrid_idx = get_specbingrid(wavegrid, 10000./sigma_in['wno'][::-1])

    sigma_array_bin = np.zeros((len(sigma_in['p']), len(sigma_in['t']), len(wavegrid)))
    for pressure_idx, pressure_val in enumerate(sigma_in['p']):
        for temperature_idx, temperature_val in enumerate(sigma_in['t']):

            print 'Compute temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)
            sigma = sigma_array[pressure_idx, temperature_idx]
            if options.sampling_method == 'geometric_average':
                # geometric average (i.e. log-average)
                sigma_rev = sigma[::-1] # reverse array
                logval = np.log(sigma_rev)
                values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(1,len(wavegrid)+1)])
                values = np.exp(values)
                values[np.isnan(values)] = 0
            elif options.sampling_method == 'algebraic_average':
                # algebraic average
                sigma_rev = sigma[::-1] # reverse array
                values = np.asarray([np.average(sigma_rev[bingrid_idx == i]) for i in xrange(1,len(wavegrid)+1)])
                values[np.isnan(values)] = 0
            elif options.sampling_method == 'random_sample':
                sigma_rev = sigma[::-1] # reverse array
                values = np.asarray([np.random.choice(sigma_rev[bingrid_idx == i], 1)[0] for i in xrange(1,len(wavegrid)+1)])
            elif options.sampling_method == 'first_sample':
                sigma_rev = sigma[::-1] # reverse array
                values = np.asarray([sigma_rev[bingrid_idx == i][0] for i in xrange(1,len(wavegrid)+1)])
            elif options.sampling_method == 'interp_sample':
                sigma_rev = sigma[::-1] # reverse array
                values = np.interp(wavegrid, 10000./sigma_in['wno'][::-1], sigma_rev)

            sigma_array_bin[pressure_idx, temperature_idx, :] = values[::-1] # reverse array back to original

    sigma_array_out = sigma_array_bin
    wno_out = 10000./wavegrid[::-1]

elif options.binning_method == 'const_wn':

    print 'Constant dwno=%.2f binning' % float(options.binning_resolution)

    # create the output wavenumber grid
    bin_wngrid = np.arange(wnmin, wnmax, float(options.binning_resolution))
    bingrid_idx = np.digitize(full_wngrid, bin_wngrid) #getting the specgrid indexes for bins

    comments.append('Constant dwno=%.2f binning' % float(options.binning_resolution))

    if options.sampling_method == 'geometric_average':
        comments.append('Linear binning: use geometric average.' )
    elif options.sampling_method == 'algebraic_average':
        comments.append('Linear binning: use algebraic average.' )
    elif options.sampling_method == 'random_sample':
        comments.append('Linear binning: use random sample.' )
    elif options.sampling_method == 'first_sample':
        comments.append('Linear binning: use first sample.' )


    sigma_array_bin = np.zeros((len(sigma_in['p']), len(sigma_in['t']), len(bin_wngrid)))
    for pressure_idx, pressure_val in enumerate(sigma_in['p']):
        for temperature_idx, temperature_val in enumerate(sigma_in['t']):

            print 'Binning: Compute temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)

            sigma = sigma_array[pressure_idx, temperature_idx]
            if options.sampling_method == 'geometric_average':
                # geometric average (i.e. log-average)
                logval = np.log(sigma)
                values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid+1))])
                values = np.exp(values)
                values[np.isnan(values)] = 0
            elif options.sampling_method == 'algebraic_average':
                # algebraic average
                values = np.asarray([np.average(sigma[bingrid_idx == i]) for i in xrange(1,len(bin_wngrid+1))])
                values[np.isnan(values)] = 0
            elif options.sampling_method == 'random_sample':
                values = np.asarray([np.random.choice(sigma[bingrid_idx == i], 1)[0] for i in xrange(1,len(bin_wngrid)+1)])
            elif options.sampling_method == 'first_sample':
                values = np.asarray([sigma[bingrid_idx == i][0] for i in xrange(1,len(bin_wngrid)+1)])
            elif options.sampling_method == 'interp_sample':
                values = np.interp(bin_wngrid, sigma_in['wno'], sigma)


            sigma_array_bin[pressure_idx, temperature_idx, :] = values
    sigma_array_out = sigma_array_bin
    wno_out = bin_wngrid
else:
    comments.append('No binning of the cross section was performed.')
    print 'No binning of the cross section was performed.'
    sigma_array_out =  sigma_array
    wno_out = full_wngrid



# interpolate to new temperature and pressure grids. Note that it is better to interpolate to the temperature
# grid in log space, and to the new pressure grid in linear space... so split interpolation in two steps
# note that the interpolation of pressures is not great, so it is better to keep a fine grid.

# interpolation of temperatures in log space
if options.temperature_list:
    sigma_array_tmp = np.zeros((len(sigma_in['p']), len(temperatures), len(wno_out)))
    comments.append('Interpolation to new temperature grid in log space')
    print 'Interpolation to new temperature grid in log space'
    for pressure_idx, pressure_val in enumerate(sigma_in['p']):
        for temperature_idx, temperature_val in enumerate(temperatures):
            print 'Interpolate temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)
            for wno_idx, wno_val in enumerate(wno_out):
                sigma_array_tmp[pressure_idx, temperature_idx, wno_idx] = np.exp(np.interp(temperature_val, sigma_in['t'], np.log(sigma_array_out[pressure_idx,:,wno_idx])))
                if np.isnan(sigma_array_tmp[pressure_idx, temperature_idx, wno_idx]):
                    sigma_array_tmp[pressure_idx, temperature_idx, wno_idx] = 0
else:
    sigma_array_tmp = sigma_array_out

# interpolation of pressures in linear space
if options.pressure_list:
    comments.append('Interpolation to new pressure grid in linear space')
    print 'Interpolation to new pressure grid in linear space'
    sigma_array = np.zeros((len(pressures), len(temperatures), len(wno_out)))
    for pressure_idx, pressure_val in enumerate(pressures):
        for temperature_idx, temperature_val in enumerate(temperatures):
            print 'Interpolate temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)
            for wno_idx, wno_val in enumerate(wno_out):
                sigma_array[pressure_idx, temperature_idx, wno_idx] = np.interp(pressure_val, sigma_in['p'], sigma_array_tmp[:,temperature_idx,wno_idx])
else:
    sigma_array = sigma_array_tmp

sigma_out = {
    'name': options.molecule_name,
    'p': pressures,
    't': temperatures,
    'wno': wno_out,
    'xsecarr': sigma_array,
    'comments': comments
}

pickle.dump(sigma_out, open(options.output_filename, 'wb'))