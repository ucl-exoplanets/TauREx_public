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
from scipy.inteprolate import interp2d
from time import gmtime, strftime




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

options, remainder = parser.parse_args()

if not options.source_sigma or \
   not options.output_filename or \
   not options.pressure_list or \
   not options.temperature_list:
   print 'Wrong input. Retry...'
   exit()

sigma_in = pickle.load(open(options.source_sigma))

# define the wavenumber range of the output cross sections
# if these limits extend above/below the input xsec, set to 0
if not options.wavenumber_range:
    wmin = np.min(sigma_in['wno'])
    wmax = np.max(sigma_in['wno'])
else:
    wnmin, wnmax = float(options.wavenumber_range.split(','))

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

# firstly reinterpolate the highres cross section to the requested pressure/temperature grid
# note that the wavenumber grid is arbitrary...
pressures = np.asarray(float(options.pressure_list.split(',')))
temperatures = np.asarray(float(options.temperature_list.split(',')))
xsecarr = sigma_in['xsecarr']
sigma_array = np.zeros((len(pressures), len(temperatures), len(sigma_in['wno'])))
for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        for wno_idx, wno_val in enumerate(temperatures):
            # to be changed with the wired inteprolation of Hill et al. 2012
            comments.append('Linear interpolation used to create the new pressure and temperature grids')
            siginterp = interp2d(sigma_in['p'], sigma_in['t'], sigma_in['xsecarr'][:,:,wno_idx])
            sigma_array[pressure_idx, temperature_idx, wno_idx] = siginterp(pressure_val, temperature_val)

# then reinterpolate to the requested wavenumber grid from wnmin to wnamx at the input resolution
full_wngrid = np.arange(wnmin, wnmax, resolution)
sigma_array_interp = np.zeros((len(pressures), len(temperatures), len(full_wngrid)))
for pressure_idx, pressure_val in enumerate(pressures):
    for temperature_idx, temperature_val in enumerate(temperatures):
        sigma_array_interp[pressure_idx, temperature_idx] =  np.interp(full_wngrid, sigma_in['wno'], sigma_array[pressure_idx, temperature_idx])

# lastly, bin if needed
if options.linear_binning:

    # create the output wavenumber grid
    bin_wngrid = np.arange(wnmin, wnmax, float(options.linear_binning))
    bingrid_idx = np.digitize(sigma_in['wno'], bin_wngrid) #getting the specgrid indexes for bins

    comments.append('Linear binning, at %f wavenumber.' % float(options.linear_binning))

    sigma_array_bin = np.zeros((len(pressures), len(temperatures), len(bin_wngrid)))
    for pressure_idx, pressure_val in enumerate(pressures):
        for temperature_idx, temperature_val in enumerate(temperatures):

            sigma = sigma_array[pressure_idx, temperature_idx]
            if options.binning_method == 'geometric_average':
                # geometric average (i.e. log-average)
                logval = np.log(sigma)
                values = np.asarray([np.average(logval[bingrid_idx == i]) for i in xrange(0,len(bin_wngrid))])
                values[np.isnan(values)] = 0
                values = np.exp(values)
            elif options.binning_method == 'algebric_average':
                # linear average
                values = np.asarray([np.average(sigma[bingrid_idx == i]) for i in xrange(0,len(bin_wngrid))])
                values[np.isnan(values)] = 0

            sigma_array[pressure_idx, temperature_idx, :] = values
    sigma_array_out = sigma_array_bin
    wno_out = bin_wngrid
else:
    comments.append('No binning of the cross section was performed.')
    sigma_array_out =  sigma_array_interp
    wno_out = full_wngrid

sigma_out = {
    'name': options.molecule_name,
    'p': pressures,
    't': temperatures,
    'wno': out_wngrid,
    'xsecarr': sigma_array_out,
    'comments': comments
}

pickle.dump(sigma_out, open(options.output_filename), 'wb')