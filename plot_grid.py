'''

Create varoius plots from the grid table

'''

import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib import rc
import optparse
import pymultinest
import os
import sys
import numpy as np
import csv
import scipy.ndimage

mpl.rcParams['axes.linewidth'] = 0.5 #set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 8})

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

parser = optparse.OptionParser()
parser.add_option('-i', '--input',
                  dest='input',
                  default=False,
)
parser.add_option('-o', '--output',
                  dest='output',
                  default='.',
)
parser.add_option('-n', '--name',
                  dest='name',
                  default=False,
)

options, remainder = parser.parse_args()

if not options.input or not options.output:
    print parser.print_help()
    exit()

# load csv table
file = open(options.input, 'rb')
reader = csv.reader(file)
grid = {}

for row in reader:

    param = row[0]
    input = float(row[1])
    res = float(row[2])
    snr = float(row[3])
    mode = int(row[4])

    value = float(row[5+mode*3])
    error = float(row[6+mode*3])
    relative_error = float(row[7+mode*3])

    if not res in grid:
        grid[res] = {}
    if not snr in grid[res]:
        grid[res][snr] = {}

    grid[res][snr][param] = (input, value, error, relative_error)

resolutions = sorted(grid.keys())
snrs = sorted(grid[grid.keys()[0]].keys())
params = grid[grid.keys()[0]][grid[grid.keys()[0]].keys()[0]].keys()


# create contour plots

res_floats = [float(i) for i in resolutions]
snr_floats = [float(i) for i in snrs]

X, Y = np.meshgrid(res_floats, snr_floats)
data_grid = np.zeros((len(resolutions), len(snrs)))
fig = plt.figure(figsize=(15,10))

p = 0 # plot number
i = 0
for param in params: # loop over all parameters

    data_grid = np.zeros((len(resolutions), len(snrs)))
    i = 0
    for res in resolutions:
        j = 0
        for snr in snrs:
            input, value, error, relative_error = grid[res][snr][param]
            data_grid[j,i] = float(relative_error)*100
            j += 1
        i += 1

    #smooth
    data_grid = scipy.ndimage.gaussian_filter(data_grid, sigma=1.0, order=0)

    p += 1 # plot number
    ax = fig.add_subplot(2,3,p)
    cont = ax.contourf(X, Y, data_grid, cmap=cm.winter)
    #ax.clabel(cont, fmt = '%2.1f', colors = 'w', fontsize=12)
    cbar = fig.colorbar(cont)
    cbar.set_label('Relative error in retrieved parameter (\%)')
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Signal-to-Noise')
    ax.set_title('%s' % param)


if options.name:
    h = fig.suptitle(options.name, fontsize=22)

plt.savefig(os.path.join(os.path.split(options.input)[0], 'contour_plots.pdf'))


#plt.savefig(os.path.join(options.folder, 'contourf_plots.pdf'))

# create surface plots

# X, Y = np.meshgrid(db['resolutions'], db['snrs'])
# data_grid= np.zeros((len(db['resolutions']), len(db['snrs'])))
# fig = plt.figure(figsize=(15,10))
# p = 0 # plot number
#
# i = 0
# for res in db['resolutions']: # loop for temperature
#     j = 0
#     for snr in db['snrs']:
#         fit = db['model_fitting'][res][snr]
#         data_grid[i,j] = float(fit['T_std']/fit['T'])*100.
#         j += 1
#     i += 1
# ax = fig.add_subplot(2,3,p+1, projection='3d')
# ax.plot_surface(X,Y, data_grid, alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
# ax.set_xlabel('Resolution')
# ax.set_ylabel('SNR')
# ax.set_zlabel(r'Normalised error bar (\%)')
# ax.set_title('Temperature')
# ax.view_init(elev=21., azim=55.)
# for mol in range(len(db['molecules_input'])): # loop over molecules
#     data_grid= np.zeros((len(db['resolutions']), len(db['snrs'])))
#     i = 0
#     for res in db['resolutions']:
#         j = 0
#         for snr in db['snrs']:
#             fit = db['model_fitting'][res][snr]
#             data_grid[i,j] = float(fit['X_std'][mol]/fit['X'][mol][0])*100.
#             j += 1
#         i += 1
#     p += 1 # plot number
#     ax = fig.add_subplot(2,3,p+1, projection='3d')
#     ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
#     ax.set_xlabel('Resolution')
#     ax.set_ylabel('SNR')
#     ax.set_zlabel(r'Normalised error bar (\%)')
#     ax.set_title('%s' % db['molecules_input'][mol])
#     ax.view_init(elev=21., azim=55.)
# if options.name:
#     h = fig.suptitle(options.name, fontsize=22)
# plt.savefig(os.path.join(options.folder, 'surface_plots.pdf'))

# create sensitivity plots

# as a function of snr
fig = plt.figure(figsize=(15,10))
cmap = cm.get_cmap('brg')
x = [float(i) for i in snrs]
i = 1
for param in params: # loop over parameters
    ax = fig.add_subplot(2,3,i)
    r = 0 # counter for plot colours
    for res in resolutions: # loop for temperature
        errors = []
        for snr in snrs:
            input, value, error, relative_error = grid[res][snr][param]
            errors.append(relative_error*100.)
        errors = np.asarray(errors)
        ax.plot(x, errors, label='R=%i' % int(res), c=cmap(r/10.))
        r += 1

    ax.set_xlabel('SNR')
    ax.set_ylabel(r'Normalised error bar (\%)')
    ax.set_title('%s ($x_\mathrm{in}$ = $ %s)' % (param, latex_float(input)))

    if i == 1:
        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(handles, labels)
        legend.draw_frame(False)
    i += 1

if options.name:
    h = fig.suptitle(options.name, fontsize=22)

plt.savefig(os.path.join(os.path.split(options.input)[0], 'params_vs_snr.pdf'))


# as a function of resolutions
fig = plt.figure(figsize=(15,10))
cmap = cm.get_cmap('brg')
x = [float(i) for i in resolutions]
i = 1
for param in params: # loop over parameters
    ax = fig.add_subplot(2,3,i)
    r = 0 # counter for plot colours
    for snr in snrs:
        errors = []
        for res in resolutions: # loop for temperature
            input, value, error, relative_error = grid[res][snr][param]
            errors.append(relative_error*100.)
        errors = np.asarray(errors)
        ax.plot(x, errors, label='SNR=%i' % int(snr), c=cmap(r/10.))
        r += 1

    ax.set_xlabel('Resolution')
    ax.set_ylabel(r'Normalised error bar (\%)')
    ax.set_title('%s ($x_\mathrm{in}$ = $ %s)' % (param, latex_float(input)))

    if i == 1:
        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(handles, labels)
        legend.draw_frame(False)
    i += 1

if options.name:
    h = fig.suptitle(options.name, fontsize=22)

plt.savefig(os.path.join(os.path.split(options.input)[0], 'params_vs_res.pdf'))

sys.exit()

# as a function of resolution

fig = plt.figure(figsize=(15,10))
cmap = cm.get_cmap('brg')
x = np.asarray(db['resolutions'])
i = 0
r = 0 # counter for plot colours
for snr in db['snrs']:
    ax = fig.add_subplot(2,3,i+1)
    errors = []
    for res in db['resolutions']: # loop for temperature
        fit = db['model_fitting'][res][snr]
        errors.append(float(fit['T_std'][0]/fit['T'][0])*100.)
    errors = np.asarray(errors)
    ax.plot(x, errors, label='snr=%.1f' % snr, c=cmap(r/10.))
    r += 1
i += 1
#ax.set_ylim((0, 20))
ax.set_xlabel('Resolution')
ax.set_ylabel(r'Normalised error bar (\%)')
ax.set_title('Temperature ($T_\mathrm{in} = %i\,$K)' % int(db['temperature_input']))
handles, labels = ax.get_legend_handles_labels()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels)
legend.draw_frame(False)
for mol in range(len(db['molecules_input'])): # loop over molecules
    r = 0 # counter for plot colours
    for snr in db['snrs']:
        ax = fig.add_subplot(2,3,i+1)
        #ax.text(0.05, 0.9, '$X_\mathrm{in} = %s$' % latex_float(db['mixing_input'][mol]),
        #        transform=ax.transAxes, fontsize=10)
        errors = []
        for res in db['resolutions']: # loop for temperature
            fit = db['model_fitting'][res][snr]
            error = float(fit['X_std'][0][mol]/fit['X'][0][mol][0])*100.
            errors.append(error)
        errors = np.asarray(errors)
        ax.plot(x, errors, c=cmap(r/5.))
        r += 1

    ax.set_xlabel('Resolution')
    ax.set_ylabel(r'Normalised error bar (\%)')
    ax.set_title('%s ($X_\mathrm{in} = %s$)' % (db['molecules_input'][mol], latex_float(db['mixing_input'][mol])))
    #ax.set_ylim((0, 50))
    i += 1
if options.name:
    h = fig.suptitle(options.name, fontsize=22)
plt.savefig(os.path.join(options.folder, 'params_vs_res.pdf'))

