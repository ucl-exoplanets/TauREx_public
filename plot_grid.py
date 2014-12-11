import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import pymultinest
import sys
import numpy as np
sys.path.append('./classes')
sys.path.append('./library')
import parameters,emission,transmission,output,fitting,tp_profile,data,preselector
from parameters import *

def get_from_grid(data, name, snr, res, value):
    for fit in data:
        if (fit['name'] == name) and (fit['snr'] == snr) and (fit['resolution'] == res):
             return fit[value]


grid_full = pickle.load(open('Output/create_grid/grid.pickle'))
nmol = [0, 1, 2, 3]
planet_name = 'warmneptune'

grid = []
for grid_element in grid_full:
    if grid_element['name'] == planet_name:
        grid.append(grid_element)

all_snr = []
all_res = []

for fit in grid:
    if not fit['snr'] in all_snr:
        all_snr.append(fit['snr'])
    if not fit['resolution'] in all_res:
        all_res.append(fit['resolution'])


# Build table

# get temperature
for res in all_res:
    for snr in all_snr:
        print 'Temperature, %i, %i, %i, %.2e, %.2e' % (get_from_grid(grid, planet_name, snr, res, 'temperature_input'),
                                                       res, snr,
                                                       get_from_grid(grid, planet_name, snr, res, 'T'),
                                                       get_from_grid(grid, planet_name, snr, res, 'T_std'))
# get molecules
for mol in nmol:
    data_grid = np.empty((len(all_res), len(all_snr)))
    i = 0
    for res in all_res:
        j = 0
        for snr in all_snr:
            print '%s, %.2e, %i, %i, %.2e, %.2e' % (get_from_grid(grid, planet_name, snr, res, 'molecules_input')[mol],
                                                    get_from_grid(grid, planet_name, snr, res, 'mixing_input')[mol],
                                                    res, snr,
                                                    get_from_grid(grid, planet_name, snr, res, 'X')[mol][0],
                                                    get_from_grid(grid, planet_name, snr, res, 'X_std')[mol])
            j += 1
        i += 1


print
print 'Plotting'
X, Y = np.meshgrid(all_res, all_snr)

# plot
print all_res, all_snr
data_grid = np.empty((len(all_res), len(all_snr)))
i = 0
for res in all_res:
    j = 0
    for snr in all_snr:
        data_grid[i,j] = get_from_grid(grid, planet_name, snr, res, 'T_std')
        j += 1
    i += 1
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)
ax.set_xlabel('Resolution')
ax.set_ylabel('SNR')
ax.set_zlabel('Error')
ax.set_title('%s - Temperature' % planet_name)
ax.view_init(elev=21., azim=55.)
plt.savefig('Output/%s_temperature.pdf' % planet_name)


for mol in nmol:
    data_grid = np.empty((len(all_res), len(all_snr)))
    i = 0
    for res in all_res:
        j = 0
        for snr in all_snr:
            data_grid[i,j] = get_from_grid(grid, planet_name, snr, res, 'X_std')[mol]
            j += 1
        i += 1


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.plot_surface(X,Y, data_grid,  alpha = 1, rstride=1, cstride=1, cmap=cm.winter, linewidth=0.5, antialiased=True)

    ax.set_xlabel('Resolution')
    ax.set_ylabel('SNR')
    ax.set_zlabel('Error')
    ax.set_title('%s - %s' % (planet_name, get_from_grid(grid, planet_name, snr, res, 'molecules_input')[mol]))
    ax.view_init(elev=21., azim=55.)
    plt.savefig('Output/%s_%s.pdf' % (planet_name, get_from_grid(grid, planet_name, snr, res, 'molecules_input')[mol]))

#plt.show()
