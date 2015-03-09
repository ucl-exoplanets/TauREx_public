###############################
#Library containing plotting routines that are otherwise too cumbersome for the main code
#
###############################

import  os, sys,string, glob
import pylab as pl
import numpy as np
from numpy import array
import scipy.ndimage as ndimage
from matplotlib.ticker import FuncFormatter,ScalarFormatter
from matplotlib import cm
from matplotlib import mlab

import logging

try:
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

scale_pow = 0

def plot_posteriors(fit_out, plot_name='plot', save2pdf=False, out_path=None, plot_contour=True):

        fig = _plot_posteriors(fit_out, plot_contour=plot_contour)

        if save2pdf:
                filename1 = os.path.join(out_path, '%s_posteriors.pdf' % plot_name)
                filename2 = os.path.join(out_path, '%s_posteriors.jpg' % plot_name)
                fig.savefig(filename1)
                fig.savefig(filename2)
                logging.info('Plot saved in %s and %s' % (filename1, filename2))


def plot_TP_profile(P, T_mean, T_sigma=None, fig=None, name=None, color='blue', save2pdf=False, out_path=None, linewidth=2.0):

    if fig is None: #accepting externally passed figure references
        fig = pl.figure()

    if name != 'downhill' and T_sigma is not None:
        pl.fill_betweenx(P*1e-5, T_mean-T_sigma, T_mean+T_sigma,alpha=0.3,color=color)
        pl.plot(T_mean-T_sigma,P*1e-5,'--',linewidth=linewidth,alpha=0.5,color=color)
        pl.plot(T_mean+T_sigma,P*1e-5,'--',linewidth=linewidth,alpha=0.5,color=color)

    pl.plot(T_mean,P*1e-5,linewidth=linewidth,color=color)
    pl.yscale('log')
    pl.xlabel('Temperature')
    pl.ylabel('Pressure (bar)')
    pl.gca().invert_yaxis()

    if save2pdf:
        if name is not None:
            filename = os.path.join(out_path, 'tp_profile_'+name+'.pdf')
        else:
            filename = os.path.join(out_path, 'tp_profile.pdf')
        fig.savefig(filename)
        logging.info('Plot saved in %s' % filename)

    return fig

def _plot_posteriors(fit_out, plot_contour=False,fontsize=30):

    try:
        params_sorted = sorted(fit_out[0]['fit_params'].iterkeys(),
                               key=lambda k: fit_out[0]['fit_params'][k]['sort_order'],
                               reverse=True)
    except:
        params_sorted = range(len(fit_out[0]))

    n_params = len(fit_out[0]['fit_params'])

    fig = pl.figure(figsize=(5*n_params, 5*n_params))

    seq = 0
    globalxlims = []

    for i in range(n_params):
        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        seq +=1

        ax = plot_1Dposterior(ax, fit_out, i, plot_contour)

        if i == 0:
            pl.ylabel("Prob. density",fontsize=fontsize)
        pl.xlabel(params_sorted[i],x=0.6, ha='right', fontsize=fontsize)
        globalxlims= ax.get_xlim()

        #scaling axis labels @todo labels can be improved but works now at least
        pl.rc('font', size=fontsize)
        SciFormatter = ScalarFormatter(useMathText=True,useOffset=True,useLocale=True)
        SciFormatter.set_powerlimits((-1, 4))
        ax.get_xaxis().set_major_formatter(SciFormatter)
        ax.get_xaxis().get_offset_text().set_x(0.9)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation(+40)
            tick.label.set_fontsize(fontsize)
        if i == 0:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_rotation(+40)
                tick.label.set_fontsize(fontsize)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_visible(False)

        #removing overlapping tick labels
        xticks = ax.xaxis.get_major_ticks()
        xticks[-1].label1.set_visible(False)

        for j in range(i):
            ax2 = pl.subplot(n_params, n_params, n_params * j + i + 1)
            seq += 1
            ax2.ticklabel_format(style='sci')
            plot_2Ddistribution(ax2, fit_out, i, j, suppressAxes=True, plot_contour=plot_contour)
            for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_rotation(+30)
                tick.label.set_fontsize(fontsize)
    pl.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.0,wspace=0.0)
    return fig

def plot_2Ddistribution(ax, fit_out, param1_idx, param2_idx, plot_contour=False, suppressAxes=False):

    # identify parameters from sorted dictionary
    try:
        params_sorted = sorted(fit_out[0]['fit_params'].iterkeys(),
                               key=lambda k: fit_out[0]['fit_params'][k]['sort_order'],
                               reverse=True)
    except:
        params_sorted = range(len(fit_out[0]))

    # combine chains from all fit_out
    allx = np.zeros(0)
    ally = np.zeros(0)
    #colors = ['Blues', 'Reds']
    for idx, solution in enumerate(fit_out):

        # plot 2D distributions. Use boxsum smoothing. Separate different fit_out with different colors. Not working
        # due to transparency issues
        # x = solution[params_sorted[param1_idx]]['trace']
        # y = solution[params_sorted[param2_idx]]['trace']
        # zd = grid_density_boxsum(np.min(x), np.min(y), np.max(x), np.max(y), 256, 256, zip(x, y))
        # axis.imshow(zd, origin='lower', cmap=cm.get_cmap(colors[idx]), alpha=0.6,
        #             extent=[np.min(x), np.max(x), np.min(y), np.max(y)])

        allx = np.concatenate((solution['fit_params'][params_sorted[param1_idx]]['trace'], allx))
        ally = np.concatenate((solution['fit_params'][params_sorted[param2_idx]]['trace'], ally))

    zd = grid_density_boxsum(np.min(allx), np.min(ally), np.max(allx), np.max(ally), 256, 256, zip(allx, ally))
    ax.imshow(zd, origin='lower', cmap=cm.get_cmap('Blues'), alpha=0.6,
                extent=[np.min(allx), np.max(allx), np.min(ally), np.max(ally)])

     #plot 1,2,3 sigma contours
    if plot_contour:
        for solution in fit_out:
            mean1 = solution['fit_params'][params_sorted[param1_idx]]['value']
            mean2 = solution['fit_params'][params_sorted[param2_idx]]['value']
            std1 = solution['fit_params'][params_sorted[param1_idx]]['std']
            std2 = solution['fit_params'][params_sorted[param2_idx]]['std']
            ax.axvline(x=mean1,linestyle='--',color='red')
            ax.axhline(y=mean2,linestyle='--',color='red')
            postcov = np.cov(solution['fit_params'][params_sorted[param1_idx]]['trace'],
                             solution['fit_params'][params_sorted[param2_idx]]['trace'])
            [eigval,eigvec] = np.linalg.eig(postcov)
            postangle = np.arctan2(eigvec[1,:],eigvec[0,:])
            eX1,eY1 = ellipse(std1,std2,postangle[0],mean1,mean2)
            eX2,eY2 = ellipse(2.0*std1,2.0*std2,postangle[0],mean1,mean2)
            eX3,eY3 = ellipse(3.0*std1,3.0*std2,postangle[0],mean1,mean2)
            ax.plot(eX1,eY1,linewidth=1.0,color='red',linestyle='-')
            ax.plot(eX2,eY2,linewidth=1.0,color='orange',linestyle='-')
            ax.plot(eX3,eY3,linewidth=1.0,color='yellow',linestyle='-')

    ax.set_xlim([np.min(allx),np.max(allx)])
    ax.set_ylim([np.min(ally),np.max(ally)])
    ax.set_aspect('auto', adjustable='box')

    if suppressAxes:
        pl.yticks([])
        pl.xticks([])

def plot_1Dposterior(ax, fit_out, param_idx, plot_contour=False):

    # identify parameters from sorted dictionary
    try:
        params_sorted = sorted(fit_out[0]['fit_params'].iterkeys(),
                               key=lambda k: fit_out[0]['fit_params'][k]['sort_order'],
                               reverse=True)
    except:
        params_sorted = range(len(fit_out[0]))
    allchains = np.zeros(0)
    for solution in fit_out:
        if plot_contour:
            mean1 = solution['fit_params'][params_sorted[param_idx]]['value']
            ax.axvline(x=mean1,linestyle='--',color='red')
        chain = solution['fit_params'][params_sorted[param_idx]]['trace']
        ax.hist(chain,
                  bins=50,normed=True,alpha=0.8, linewidth=2.0,
                  histtype='step', color='#408DBD') # @todo custom bins?
        allchains = np.concatenate((allchains, chain))
        ax.set_xlim((np.min(allchains), np.max(allchains)))
    return ax

def exp_formatter_fun(x, p):
    global scale_pow
    if x < 1.0:
        return "%.2f" % (x * (10 ** np.abs(scale_pow)))
    if x >= 1.0:
        return "%.2f" % (x)

def ellipse(ra,rb,ang,x0,y0,Nb=100):
    # Define a function to make the ellipses
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=np.linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y

def grid_density_boxsum(x0, y0, x1, y1, w, h, data):
    kx = (w - 1) / (x1 - x0)
    ky = (h - 1) / (y1 - y0)
    r = 15
    border = r * 2
    imgw = (w + 2 * border)
    imgh = (h + 2 * border)
    img = [0] * (imgw * imgh)
    for x, y in data:
        ix = int((x - x0) * kx) + border
        iy = int((y - y0) * ky) + border
        if 0 <= ix < imgw and 0 <= iy < imgh:
            img[iy * imgw + ix] += 1
    for p in xrange(4):
        boxsum(img, imgw, imgh, r)
    a = np.array(img).reshape(imgh,imgw)
    b = a[border:(border+h),border:(border+w)]
    return b

def boxsum(img, w, h, r):
    st = [0] * (w+1) * (h+1)
    for x in xrange(w):
        st[x+1] = st[x] + img[x]
    for y in xrange(h):
        st[(y+1)*(w+1)] = st[y*(w+1)] + img[y*w]
        for x in xrange(w):
            st[(y+1)*(w+1)+(x+1)] = st[(y+1)*(w+1)+x] + st[y*(w+1)+(x+1)] - st[y*(w+1)+x] + img[y*w+x]
    for y in xrange(h):
        y0 = max(0, y - r)
        y1 = min(h, y + r + 1)
        for x in xrange(w):
            x0 = max(0, x - r)
            x1 = min(w, x + r + 1)
            img[y*w+x] = st[y0*(w+1)+x0] + st[y1*(w+1)+x1] - st[y1*(w+1)+x0] - st[y0*(w+1)+x1]





# def _plot_mcmc_results(DATA, parameters=False, save2pdf=False, out_path=None):
#     #routine plotting MCMC posterior distributions and marginals
#     # n_params = len(parameters)
#
#     names = DATA.stats().keys()
#     n_params = len(names)
#
#     if not parameters:
#         parameters = names
#     else:
#         parameters = parameters
#
#     fig = pl.figure(figsize=(5*n_params, 5*n_params))
#     pl.title('MCMC posteriors')
#     for i in range(n_params):
#         ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
#         pl.hist(DATA.trace(names[i])[:],color='k')
#         pl.ylabel("Probability")
#         pl.xlabel(parameters[i])
#         for tick in ax.xaxis.get_major_ticks():
#             tick.label.set_rotation(-30)
#         for j in range(i):
#             ax = pl.subplot(n_params, n_params, n_params * j + i + 1)
#             ax.ticklabel_format(style='sci')
#             pl.plot(DATA.trace(names[i])[:],DATA.trace(names[j])[:],'.',color='k')
#             pl.xlabel(parameters[i])
#             pl.ylabel(parameters[j])
#             for tick in ax.xaxis.get_major_ticks():
#                     tick.label.set_rotation(-30)
#
#     if save2pdf:
#             filename = os.path.join(out_path, 'mcmc_posteriors.pdf')
#             fig.savefig(filename)
#             logging.info('Plot saved in %s' % filename)
#
# def plot_mcmc_results(mcmcdir, parameters,  save2pdf=False, out_path=None, plot_contour=True):
#     #routine plotting MCMC chains. @todo NOT WORKING
#
#     n_params = len(parameters) #number of parameters
#
#     mcmc = {}
#     for thread in glob.glob(os.path.join(mcmcdir, '*')):
#         print thread
#         print thread[-1]
#         print
#
#         id = int(thread[-1])
#         mcmc[id] = {}
#
#         with open(os.path.join(mcmcdir, 'thread_0', 'state.txt'),'r') as statefile:
#             mcmc[id]['state'] = eval(statefile.read())
#
#         chainlist = glob.glob(os.path.join(thread, 'Chain_0', '*'))
#         for trace in chainlist:
#             traceid = string.rsplit(trace, '/')[-1][:-4]
#             mcmc[id][traceid] ={}
#             mcmc[id][traceid]['data'] = np.loadtxt(trace)
#             mcmc[id][traceid]['stats'] = {}
#             mcmc[id][traceid]['stats'][0] = {}
#             mcmc[id][traceid]['stats'][0]['sigma'] = np.std(mcmc[id][traceid]['data'])
#             mcmc[id][traceid]['stats'][0]['mean'] = np.mean(mcmc[id][traceid]['data'])
#
#     #     print mcmc.keys()
#     params_mcmc = mcmc[0]['state']['stochastics'].keys()
#     params_mcmc[1:] = np.sort(params_mcmc[1:])
#
#     fig_mcmc = plot_posteriors(mcmc,n_params,params_mcmc,parameters,plot_contour=plot_contour,alpha=0.05)
#
#     if save2pdf:
#             filename = os.path.join(out_path, 'mcmc_posteriors.pdf')
#             fig_mcmc.savefig(filename)
#             logging.info('Plot saved in %s' % filename)
#
#
