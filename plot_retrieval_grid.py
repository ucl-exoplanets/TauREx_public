'''

Create varoius plots from the grid table

'''

import pickle
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import cm

#from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib import rc
import optparse
import os
import sys
import numpy as np
import csv
import scipy.ndimage

sys.path.append('./classes')
import parameters
from parameters import *


mpl.rcParams['axes.linewidth'] = 0.5 #set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 10})

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  help='Parameter filename')
parser.add_option('-l', '--levels',
                  dest='levels',
                  default=False
)
parser.add_option('-b', '--black',
                  dest='blackbg',
                  default=None
)
options, remainder = parser.parse_args()

if options.levels:
    levels = [int(m) for m in options.levels.split(',')]
else:
    levels = [0, 10, 20, 30, 40, 50, 60, 70]

cmaps = ['gist_earth'] #, 'terrain', 'ocean', 'gist_stern', 'brg', 'cubehelix', 'autumn', 'summer', 'winter']
params_sorted = ['1H2-16O', '12C-1H4', '12C-16O', '14N-1H3']
xin = ['10^{-4}', '10^{-4}', '10^{-4}', r'0.5 \times 10^{-4}']

#initialising parameters object
params = parameters(options.param_filename)

grid = {}
resolutions = []
errors = []
for res in os.listdir(params.out_path):
    if os.path.isdir(os.path.join(params.out_path, res)):
        if not int(res[1:]) in grid:
            grid[int(res[1:])] = {}
            resolutions.append(int(res[1:]))
        for error in os.listdir(os.path.join(params.out_path, res)):
            if os.path.isdir(os.path.join(params.out_path, res, error)):
                if not error in grid[int(res[1:])]:
                    grid[int(res[1:])][int(error[1:])] = {}
                    errors.append(int(error[1:]))
                filename = os.path.join(params.out_path, res, error, 'NEST_out.db')
                if os.path.isfile(filename):
                    out = pickle.load(open(filename, 'rb'))
                    for param in out[0]['fit_params']:
                        value = np.power(10, out[0]['fit_params'][param]['value'])
                        value_error = value * np.log(10) * out[0]['fit_params'][param]['std']
                        relative_error = value_error/value
                        grid[int(res[1:])][int(error[1:])][param] = (value, value_error, relative_error)

resolutions = list(set(resolutions))
errors = list(set(errors))

params = out[0]['fit_params']

X, Y = np.meshgrid(resolutions, errors)
data_grid = np.zeros((len(resolutions), len(errors)))

for cmap in cmaps:

    # create contour plots

    fig = plt.figure(1, figsize=(5,4))

    plot_grid = AxesGrid(fig, 111,
                    nrows_ncols = (2, 2),
                    aspect=False,
                    axes_pad = 0.3,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="5%",
                    cbar_pad=0.2)


    p = 0 # plot number
    i = 0

    for param in params_sorted: # loop over all parameters

        data_grid = np.zeros((len(resolutions), len(errors)))

        i = 0
        for res in resolutions:
            j = 0
            for error in errors:
                value, error, relative_error = grid[res][error][param]
                data_grid[j,i] = float(relative_error)*100
                j += 1
            i += 1

        #smooth
        data_grid = scipy.ndimage.gaussian_filter(data_grid, sigma=1.0, order=0)

        #im = plot_grid[p].contourf(X, Y, data_grid, levels=levels, cmap=cmap)
        im = plot_grid[p].contourf(X, Y, data_grid, cmap=cmap)
        plot_grid[p].set_xlabel('Resolution')
        plot_grid[p].set_ylabel('Error (ppm)')
        if param == '14N-1H3':
            title = 'NH$_3$'
        if param == '1H2-16O':
            title = 'H$_2$O'
        if param == '12C-16O':
            title = 'CO'
        if param == '12C-1H4':
            title = 'CH$_4$'

        print data_grid

        plot_grid[p].set_title(r'%s ($X_\mathrm{in} = %s$)' % (title, xin[p]), fontsize=10)

        p += 1 # plot number


    cbar = plt.colorbar(im, cax=plot_grid.cbar_axes[0])
    cbar.ax.set_ylabel('Relative error in retrieved parameter (\%)')

    plot_grid.axes_llc.set_xticks([40, 60, 80])
    plot_grid.axes_llc.set_yticks([50, 100, 150])

    plt.tight_layout()
    #plt.savefig(os.path.join(output, '%scontour_%s.pdf' % (name, cmap)))
    plt.show()
    plt.clf()
