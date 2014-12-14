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

mpl.rcParams['axes.linewidth'] = 0.5 #set the value globally
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 8})


sys.path.append('./classes')
sys.path.append('./library')
import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

parser = optparse.OptionParser()
parser.add_option('-f', '--folder',
                  dest='folder',
                  default=False,
)
parser.add_option('-n', '--name',
                  dest='name',
                  default=False,
)

options, remainder = parser.parse_args()

if not options.folder:
    print parser.print_help()
    exit()

# load database
db = pickle.load(open(os.path.join(options.folder, 'grid.db'), 'rb'))

# Build table
grid_filename = os.path.join(options.folder, 'grid.csv')


f = open(grid_filename,'w')
for res in db['resolutions']: # loop over molecules
    for snr in db['snrs']:
        fit = db['model_fitting'][res][snr]
        f.write('Temperature, %i, %i, %i, %.2e, %.2e \n' % (db['temperature_input'],res, snr, fit['T'], fit['T_std']))
for mol in range(len(db['molecules_input'])): # loop for each molecule
    for res in db['resolutions']:
        for snr in db['snrs']:
            fit = db['model_fitting'][res][snr]
            f.write('%s, %.2e, %i, %i, %.2e, %.2e \n' % (db['molecules_input'][mol],
                                                         db['mixing_input'][mol],
                                                         res, snr,
                                                         fit['X'][mol][0],
                                                         fit['X_std'][mol]))
f.close()

# create surface plots
X, Y = np.meshgrid(db['resolutions'], db['snrs'])
data_grid= np.empty((len(db['resolutions']), len(db['snrs'])))

fig = plt.figure(figsize=(15,10))
p = 0 # plot number
i = 0
for res in db['resolutions']: # loop for temperature
    j = 0
    for snr in db['snrs']:
        fit = db['model_fitting'][res][snr]
        data_grid[i,j] = float(fit['T_std']/fit['T'])*100.
        j += 1
    i += 1
i = 0
fig = plt.figure()
ax = fig.add_subplot(2,3,p+1, projection='3d')
p += 1
ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
ax.set_xlabel('Resolution')
ax.set_ylabel('SNR')
ax.set_zlabel(r'Normalised error bar (\%)')
ax.set_title('Temperature')
ax.view_init(elev=21., azim=55.)

for mol in range(len(db['molecules_input'])): # loop over molecules
    data_grid= np.empty((len(db['resolutions']), len(db['snrs'])))
    i = 0
    for res in db['resolutions']:
        j = 0
        for snr in db['snrs']:
            fit = db['model_fitting'][res][snr]
            data_grid[i,j] = float(fit['X_std'][mol]/fit['X'][mol][0])*100.
        j += 1
    i += 1
    fig = plt.figure()
    ax = fig.add_subplot(2,3,p+1, projection='3d')
    p += 1 # plot number
    ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('SNR')
    ax.set_zlabel(r'Normalised error bar (\%)')
    ax.set_title('%s' % db['molecules_input'][mol])
    ax.view_init(elev=21., azim=55.)
if options.name:
    h = fig.suptitle(options.name, fontsize=22)
plt.savefig(os.path.join(options.folder, 'surface_plots.pdf'))

# create sensitivity plots
fig = plt.figure(figsize=(15,10))
cmap = cm.get_cmap('winter')
i = 0
r = 0 # counter for plot colours
for res in db['resolutions']: # loop for temperature
    ax = fig.add_subplot(2,3,i+1)
    ax.text(0.05, 0.9, '$T_\mathrm{in} = %i\,$K' % int(db['temperature_input']),
            transform=ax.transAxes, fontsize=10)
    errors = []
    for snr in db['snrs']:
        fit = db['model_fitting'][res][snr]
        errors.append(float(fit['T_std']/fit['T'])*100.)
    ax.plot(db['snrs'], errors, label='R=%i' % int(res), c=cmap(r/5.))
    r += 1
i += 1
ax.set_ylim((0, 50))
ax.set_xlabel('SNR')
ax.set_ylabel(r'Normalised error bar (\%)')
ax.set_title('Temperature')
handles, labels = ax.get_legend_handles_labels()
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels)
legend.draw_frame(False)

for mol in range(len(db['molecules_input'])): # loop over molecules
    r = 0 # counter for plot colours
    for res in db['resolutions']: # loop for temperature
        ax = fig.add_subplot(2,3,i+1)
        ax.text(0.05, 0.9, '$X_\mathrm{in} = %s$' % latex_float(db['mixing_input'][mol]),
                transform=ax.transAxes, fontsize=10)

        errors = []
        for snr in db['snrs']:
            fit = db['model_fitting'][res][snr]
            error = float(fit['X_std'][mol]/fit['X'][mol][0])*100.
            errors.append(error)
        ax.plot(db['snrs'], errors, c=cmap(r/5.))
        r += 1

    ax.set_xlabel('SNR')
    ax.set_ylabel(r'Normalised error bar (\%)')
    ax.set_title('%s' % db['molecules_input'][mol])
    ax.set_ylim((0, 200))
    i += 1

if options.name:
    h = fig.suptitle(options.name, fontsize=22)

plt.savefig(os.path.join(options.folder, 'params_vs_snr.pdf'))
