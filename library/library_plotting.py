###############################
#Library containing plotting routines that are otherwise too cumbersome for the main code
#
###############################

import  os, sys,string
import pylab as pl
import numpy as np
from numpy import array
import scipy.ndimage as ndimage
from matplotlib.ticker import FuncFormatter
import logging

try:  
    import pymultinest
    multinest_import = True
except:
    multinest_import = False


def plot_multinest_results(DATA, parameters,  save2pdf=False, out_path=None):

    #routine plotting multinest posterior distributions and marginals
    n_params = len(parameters)
#         s = a.get_stats()
    p = pymultinest.PlotMarginalModes(DATA)
    fig = pl.figure(figsize=(5*n_params, 5*n_params))
    pl.title('Multinest posteriors')
    #plt.subplots_adjust(wspace=0, hspace=0)
#         fig.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    for i in range(n_params):

        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        p.plot_marginal(i, with_ellipses=True, with_points=False, grid_points=50)

#             p1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#             plt.xticks(rotation=70)

        pl.ylabel("Probability")
        pl.xlabel(parameters[i])
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation(-30)

        for j in range(i):
            ax = pl.subplot(n_params, n_params, n_params * j + i + 1)
            ax.ticklabel_format(style='sci')
            #plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
            p.plot_conditional(i, j, with_ellipses=False, with_points=True, grid_points=50)
            pl.xlabel(parameters[i])
            pl.ylabel(parameters[j])
            for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_rotation(-30)
    if save2pdf:
            filename = os.path.join(out_path, 'nested_posteriors.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)


def plot_mcmc_results(DATA, parameters=False, save2pdf=False, out_path=None):
    #routine plotting MCMC posterior distributions and marginals
    # n_params = len(parameters)

    names = DATA.stats().keys()
    n_params = len(names)

    if not parameters:
        parameters = names
    else:
        parameters = parameters

    fig = pl.figure(figsize=(5*n_params, 5*n_params))
    pl.title('MCMC posteriors')
    for i in range(n_params):
        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        pl.hist(DATA.trace(names[i])[:],color='k')
        pl.ylabel("Probability")
        pl.xlabel(parameters[i])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation(-30)
        for j in range(i):
            ax = pl.subplot(n_params, n_params, n_params * j + i + 1)
            ax.ticklabel_format(style='sci')
            pl.plot(DATA.trace(names[i])[:],DATA.trace(names[j])[:],'.',color='k')
            pl.xlabel(plotnames[i])
            pl.ylabel(plotnames[j])
            for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_rotation(-30)

    if save2pdf:
            filename = os.path.join(self.params.out_path, 'mcmc_posteriors.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)


