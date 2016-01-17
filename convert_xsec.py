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

import sys, os, optparse, string, pickle, glob
import numpy as np
#from scipy.inteprolate import interp2d
from time import gmtime, strftime

import matplotlib.pylab as plt


parser = optparse.OptionParser()
parser.add_option('-s', '--source_sigma',
                  dest='source_sigma',
                  default=None,
)
parser.add_option('-o', '--output_filename',
                  dest='output_filename',
                  default='sigma_array.db',
)
parser.add_option('-p', '--pressure_list',
                  dest='pressure_list',
                  default=None,
                  )
parser.add_option('-t', '--temperature_list',
                  dest='temperature_list',
                  default=None,
                  )
parser.add_option('-w', '--wavenumber_range',
                  dest='wavenumber_range',
                  default=None,
                  )
parser.add_option('-b', '--binning_resolution',
                  dest='linear_binning',
                  default=None,
)
parser.add_option('-m', '--binning_method',
                  dest='binning_method',
                  default='geometric_average',
)
parser.add_option('-r', '--input_resolution',
                  dest='input_resolution',
                  default=None,
)
parser.add_option('-n', '--molecule_name',
                  dest='molecule_name',
                  default=None,
)

options, remainder = parser.parse_args()

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
        sigma_array[pressure_idx, temperature_idx] =  np.interp(full_wngrid, sigma_in['wno'], sigma_in['xsecarr'][pressure_idx, temperature_idx])

# interpolate to new temperature and pressure grids. Note that it is better to interpolate to the temperature
# grid in log space, and to the new pressure grid in linear space... so split interpolation in two steps
# note that the interpolation of pressures is not great, so it is better to keep a fine grid.

# interpolation of temperatures in log space
sigma_array_tmp = np.zeros((len(sigma_in['p']), len(temperatures), len(full_wngrid)))
comments.append('Interpolation to new temperature grid in log space')
print 'Interpolation to new temperature grid in log space'
for pressure_idx, pressure_val in enumerate(sigma_in['p']):
    for temperature_idx, temperature_val in enumerate(temperatures):
        print 'Interpolate temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)
        for wno_idx, wno_val in enumerate(full_wngrid):
            sigma_array_tmp[pressure_idx, temperature_idx, wno_idx] = np.exp(np.interp(temperature_val, sigma_in['t'], np.log(sigma_array[pressure_idx,:,wno_idx])))
            if np.isnan(sigma_array_tmp[pressure_idx, temperature_idx, wno_idx]):
                sigma_array_tmp[pressure_idx, temperature_idx, wno_idx] = 0

# interpolation of pressures in linear space
comments.append('Interpolation to new pressure grid in linear space')
sigma_array = np.zeros((len(pressures), len(temperatures), len(full_wngrid)))
for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        print 'Interpolate temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)
        for wno_idx, wno_val in enumerate(full_wngrid):
            sigma_array[pressure_idx, temperature_idx, wno_idx] = np.interp(pressure_val, sigma_in['p'], sigma_array_tmp[:,temperature_idx,wno_idx])

# lastly, bin if needed.
# For now, only linear binning available. Ideally implement optimal binning for input resolution.
if options.linear_binning:

    print 'Linear binning'

    # create the output wavenumber grid
    bin_wngrid = np.arange(wnmin, wnmax, float(options.linear_binning))
    bingrid_idx = np.digitize(full_wngrid, bin_wngrid) #getting the specgrid indexes for bins

    comments.append('Linear binning, at %f wavenumber resolution.' % float(options.linear_binning))

    if options.binning_method == 'geometric_average':
        comments.append('Linear binning: use geometric average.' )
    elif options.binning_method == 'algebraic_average':
        comments.append('Linear binning: use algebraic average.' )

    sigma_array_bin = np.zeros((len(pressures), len(temperatures), len(bin_wngrid)))
    for pressure_idx, pressure_val in enumerate(pressures):
        for temperature_idx, temperature_val in enumerate(temperatures):

            print 'Compute temperature %.1f, pressure %.3e' % (temperature_val, pressure_val)

            sigma = sigma_array[pressure_idx, temperature_idx]
            if options.binning_method == 'geometric_average':
                # geometric average (i.e. log-average)
                logval = np.log(sigma)
                values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(0,len(bin_wngrid))])
                values = np.exp(values)
                values[np.isnan(values)] = 0
            elif options.binning_method == 'algebraic_average':
                # algebraic average
                values = np.asarray([np.average(sigma[bingrid_idx == i]) for i in xrange(0,len(bin_wngrid))])
                values[np.isnan(values)] = 0
            # here you can try other methods, such as random sampling


            sigma_array_bin[pressure_idx, temperature_idx, :] = values
    sigma_array_out = sigma_array_bin
    wno_out = bin_wngrid
else:
    comments.append('No binning of the cross section was performed.')
    sigma_array_out =  sigma_array
    wno_out = full_wngrid

sigma_out = {
    'name': options.molecule_name,
    'p': pressures,
    't': temperatures,
    'wno': wno_out,
    'xsecarr': sigma_array_out,
    'comments': comments
}

pickle.dump(sigma_out, open(options.output_filename, 'wb'))