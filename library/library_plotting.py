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

def plot_posteriors(fit_out, params_names=[], plot_name='plot', save2pdf=False, out_path=None, plot_contour=True, color='Blues'):

        fig = _plot_posteriors(fit_out, params_names=params_names, plot_contour=plot_contour, color=color)

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

def _plot_posteriors(fit_out, params_names=[], plot_contour=False,fontsize=30, color='Blues'):

    if params_names == []:
        params_names = fit_out[0]['fit_params'].keys()

    n_params = len(params_names)

    fig = pl.figure(figsize=(5*n_params, 5*n_params))

    seq = 0
    globalxlims = []

    for i, param_name in enumerate(params_names):

        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        seq +=1

        ax = plot_1Dposterior(ax, fit_out, params_names, i, plot_contour)

        if i == 0:
            pl.ylabel("Prob. density",fontsize=fontsize)
        pl.xlabel(param_name, x=0.6, ha='right', fontsize=fontsize)
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
            plot_2Ddistribution(ax2, fit_out, params_names, i, j, suppressAxes=True, plot_contour=plot_contour, color=color)
            for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_rotation(+30)
                tick.label.set_fontsize(fontsize)
    pl.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.0,wspace=0.0)
    return fig

def plot_2Ddistribution(ax, fit_out, params_names, param1_idx, param2_idx, plot_contour=False, suppressAxes=False, color='Blues'):


    # combine chains from all fit_out
    allx = np.zeros(0)
    ally = np.zeros(0)

    #colors = ['Blues', 'Reds']
    for idx, solution in enumerate(fit_out):

        # plot 2D distributions. Use boxsum smoothing. Separate different fit_out with different colors. ** Not working
        # due to transparency issues
        # x = solution[params_names[param1_idx]]['trace']
        # y = solution[params_names[param2_idx]]['trace']
        # zd = grid_density_boxsum(np.min(x), np.min(y), np.max(x), np.max(y), 256, 256, zip(x, y))
        # axis.imshow(zd, origin='lower', cmap=cm.get_cmap(colors[idx]), alpha=0.6,
        #             extent=[np.min(x), np.max(x), np.min(y), np.max(y)])

        # concatenate all separate traces to get a single 2D distribution
        allx = np.concatenate((solution['fit_params'][params_names[param1_idx]]['trace'], allx))
        ally = np.concatenate((solution['fit_params'][params_names[param2_idx]]['trace'], ally))

    zd = grid_density_boxsum(np.min(allx), np.min(ally), np.max(allx), np.max(ally), 256, 256, zip(allx, ally))
    ax.imshow(zd, origin='lower', cmap=cm.get_cmap(color),
              extent=[np.min(allx), np.max(allx), np.min(ally), np.max(ally)])

     #plot 1,2,3 sigma contours
    if plot_contour:
        for solution in fit_out:
            mean1 = solution['fit_params'][params_names[param1_idx]]['value']
            mean2 = solution['fit_params'][params_names[param2_idx]]['value']
            std1 = solution['fit_params'][params_names[param1_idx]]['std']
            std2 = solution['fit_params'][params_names[param2_idx]]['std']
            ax.axvline(x=mean1,linestyle='--',color='red')
            ax.axhline(y=mean2,linestyle='--',color='red')
            if std1 > 0 and std2 > 0:
                postcov = np.cov(solution['fit_params'][params_names[param1_idx]]['trace'],
                                 solution['fit_params'][params_names[param2_idx]]['trace'])
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

def plot_1Dposterior(ax, fit_out, params_names, param_idx, plot_contour=False):

    allchains = np.zeros(0)
    for solution in fit_out:
        if plot_contour:
            mean1 = solution['fit_params'][params_names[param_idx]]['value']
            ax.axvline(x=mean1,linestyle='--',color='red')
        chain = solution['fit_params'][params_names[param_idx]]['trace']
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

