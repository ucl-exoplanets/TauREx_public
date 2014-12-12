import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import optparse
import pymultinest
import os
import sys
import numpy as np
sys.path.append('./classes')
sys.path.append('./library')
import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *

parser = optparse.OptionParser()
parser.add_option('-f', '--folder',
                  dest='folder',
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
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
ax.set_xlabel('Resolution')
ax.set_ylabel('SNR')
ax.set_zlabel('Normalised error bar (%)')
ax.set_title('Temperature')
ax.view_init(elev=21., azim=55.)
plt.savefig(os.path.join(options.folder, 'temperature.pdf'))
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
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('SNR')
    ax.set_zlabel('Normalised error bar (%)')
    ax.set_title('%s' % db['molecules_input'][mol])
    ax.view_init(elev=21., azim=55.)
    plt.savefig(os.path.join(options.folder, '%s.pdf' % db['molecules_input'][mol]))
plt.show()