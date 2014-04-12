###############################
#Library containing plotting routines that are otherwise too cumbersome for the main code
#
###############################

import pylab
from pylab import *
import numpy 
from numpy import *
import pymultinest


def plot_multinest_results(DATA, parameters, save2pdf=False):
    #routine plotting multinest posterior distributions and marginals
    n_params = len(parameters)
#         s = a.get_stats()
    p = pymultinest.PlotMarginalModes(DATA)
    figure(figsize=(5*n_params, 5*n_params))
    title('Multinest posteriors')
    #plt.subplots_adjust(wspace=0, hspace=0)
#         fig.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    for i in range(n_params):
        ax = subplot(n_params, n_params, n_params * i + i + 1)
        p.plot_marginal(i, with_ellipses=True, with_points=False, grid_points=50)

#             p1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#             plt.xticks(rotation=70)

        ylabel("Probability")
        xlabel(parameters[i])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(-30)

    for j in range(i):
        ax = plt.subplot(n_params, n_params, n_params * j + i + 1)
        ax.ticklabel_format(style='sci')
        #plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
        p.plot_conditional(i, j, with_ellipses=False, with_points=True, grid_points=50)
        xlabel(parameters[i])
        ylabel(parameters[j])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(-30)


def plot_mcmc_results(DATA,parameters=False, save2pdf=False):
    #routine plotting MCMC posterior distributions and marginals
    # n_params = len(parameters)

    names = DATA.stats().keys()
    n_params = len(names)

    if not parameters:
        plotnames = names
    else:
        plotnames = parameters

    figure(figsize=(5*n_params, 5*n_params))
    title('MCMC posteriors')
    for i in range(n_params):
        ax = subplot(n_params, n_params, n_params * i + i + 1)
        hist(DATA.trace(names[i])[:],color='k')
        ylabel("Probability")
        xlabel(plotnames[i])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(-30)
    for j in range(i):
        ax = plt.subplot(n_params, n_params, n_params * j + i + 1)
        ax.ticklabel_format(style='sci')
        plot(DATA.trace(names[i])[:],DATA.trace(names[j])[:],'.',color='k')
        xlabel(plotnames[i])
        ylabel(plotnames[j])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(-30)


