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
import matplotlib.colors as mplot_colors
from matplotlib import mlab
import matplotlib.patches
import itertools
import scipy

import logging

try:
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

scale_pow = 0

def plot_posteriors(fit_out, params_names=[], plot_name='plot', save2pdf=False, out_path=None, plot_contour=True, color='BuPu',log_cmap=False):
        
        fig = _plot_posteriors(fit_out, plot_name=plot_name, params_names=params_names, plot_contour=plot_contour, color=color,log_cmap=log_cmap)

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

def _plot_posteriors(fit_out, plot_name=None, params_names=None, plot_contour=False,fontsize=30, color='Blues',log_cmap=False):

    if params_names == None:
        params_names = fit_out[0]['fit_params'].keys()

    n_params = len(params_names)

    fig = pl.figure(figsize=(5*n_params, 5*n_params))

    seq = 0
    globalxlims = []

    for i, param_name in enumerate(params_names):

        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        seq +=1

        if plot_name == 'NEST' or plot_name == 'NEST_clrinv':
            ax = NEST_plot_conditional(ax, fit_out, params_names, i, suppressAxes=False, plot_contour=plot_contour, color=color,log_cmap=log_cmap)
        else:
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
            if plot_name == 'NEST' or plot_name == 'NEST_clrinv':
                NEST_plot_conditional(ax2, fit_out, params_names, i, j, suppressAxes=True, plot_contour=plot_contour, color=color,log_cmap=log_cmap)
            else:
                plot_2Ddistribution(ax2, fit_out, params_names, i, j, suppressAxes=True, plot_contour=plot_contour, color=color,log_cmap=log_cmap)
            for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_rotation(+30)
                tick.label.set_fontsize(fontsize)
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.0,wspace=0.0)
    return fig


def NEST_plot_conditional(ax, fit_out, params_names, param1_idx, param2_idx=None, plot_contour=False, suppressAxes=False,  color='Blues',log_cmap=False):

    with_ellipses = True
    with_points = False
    only_interpolate = False
    use_log_values = False
    grid_points = 100
    marginalization_type='sum' # 'max', 'mean'

    # this code is largely taken from PyMultinest plot.py PlotMarginalModes.plot_conditional()

    # combine chains from all fit_out
    dim1_column = np.zeros(0)
    dim2_column = np.zeros(0)
    values = np.zeros(0)

    for idx, solution in enumerate(fit_out):
        dim1_column = np.concatenate((solution['fit_params'][params_names[param1_idx]]['trace'], dim1_column))
        if param2_idx is not None:
            dim2_column = np.concatenate((solution['fit_params'][params_names[param2_idx]]['trace'], dim2_column))
        values = np.concatenate((solution['weights'], values))

    min1 = np.min(dim1_column)
    max1 = np.max(dim1_column)
    mean1 = solution['fit_params'][params_names[param1_idx]]['value']
    std1 = solution['fit_params'][params_names[param1_idx]]['std']

    if param2_idx is not None:
        coords = np.array([dim1_column, dim2_column]).transpose()
        min2 = np.min(dim2_column)
        max2 = np.max(dim2_column)
        mean2 = solution['fit_params'][params_names[param2_idx]]['value']
        std2 = solution['fit_params'][params_names[param2_idx]]['std']
    else:
        coords = dim1_column.transpose()

    if use_log_values:
        values = np.log(values)


    # determining min/max
    # min1 = min([solution['fit_params'][params_names[param1_idx]]['value'] -
    #             3*solution['fit_params'][params_names[param1_idx]]['std'] for solution in fit_out])
    # max1 = max([solution['fit_params'][params_names[param1_idx]]['value'] +
    #             3*solution['fit_params'][params_names[param1_idx]]['std'] for solution in fit_out])
    # min2 = min([solution['fit_params'][params_names[param2_idx]]['value'] -
    #             3*solution['fit_params'][params_names[param2_idx]]['std'] for solution in fit_out])
    # max2 = max([solution['fit_params'][params_names[param2_idx]]['value'] +
    #             3*solution['fit_params'][params_names[param2_idx]]['std'] for solution in fit_out])


    n = grid_points
    binsize1 = (max1 - min1) / n

    if param2_idx is not None:
        m = n
        grid_x, grid_y = np.mgrid[min1:max1:n*1j, min2:max2:n*1j]
        binsize2 = (max2 - min2) / m
    else:
        m = 1
        grid_x = np.mgrid[min1:max1:n*1j]
        grid_y = [0]

    grid_z = np.zeros((n,m))
    minvalue = values.min()
    maxvalue = values.max()

    # for each grid item, find the matching points and put them in.
    for row, col in itertools.product(list(range(len(grid_x))), list(range(len(grid_y)))):

        if param2_idx is not None:
            xc = grid_x[row,col]
            here_x = np.abs(dim1_column - xc) < binsize1 / 2.
            yc = grid_y[row,col]
            here_y = np.abs(dim2_column - yc) < binsize2 / 2.
        else:
            xc = grid_x[row]
            here_x = np.abs(dim1_column - xc) < binsize1 / 2.
            here_y = True

        bin = values[np.logical_and(here_x, here_y)]
        if bin.size != 0:
            if marginalization_type == 'max':
                grid_z[row,col] = bin.max()
            elif marginalization_type == 'sum':
                grid_z[row,col] = bin.sum()
            elif marginalization_type == 'mean':
                grid_z[row,col] = bin.mean()
            elif marginalization_type == 'count':
                grid_z[row,col] = bin.size
            else:
                assert False, "marginalization_type should be mean, sum or max"
        else:
            grid_z[row,col] = minvalue

    if only_interpolate:
    #   version A: interpolated -- may look weird because of the
    #              loss of dimensions
        grid_z = scipy.interpolate.griddata(coords, values, (grid_x, grid_y), method='cubic')

    ax.set_xlim((min1,max1))
    if param2_idx is not None:

        #customising the colormap
        new_cmap = pl.cm.get_cmap(color)
        new_cmap.set_under('w')
        if log_cmap is True:
            cmap_norm = mplot_colors.LogNorm()
        else:
            cmap_norm = mplot_colors.Normalize()

        ax.set_ylim((min2,max2))
        ax.imshow(grid_z.transpose(), origin='lower', aspect='auto',
                 cmap=new_cmap,vmin=1e-10, norm=cmap_norm,extent=(min1,max1,min2,max2))
        #plt.colorbar()
    else:
		ax.plot(grid_x, grid_z[:,0], '-', color='grey', drawstyle='steps')

    if use_log_values:
        levels = [maxvalue, maxvalue - .5, maxvalue - 1.0, maxvalue - 2.0]
    else:
        levels = [maxvalue, maxvalue / 3, maxvalue / 10, maxvalue / 100]
    leveltitles = ['max', 'max/3', 'max/10', 'max/100']
    #ax.contour(grid_x, grid_y, grid_z, levels, linewidths=0.5, colors='k')

    if with_points and param2_idx is not None:
	    ax.scatter(dim1_column, dim2_column, marker='+', color='black', s=1, alpha=0.3)



    if plot_contour:
        for solution in fit_out:
            mean1 = solution['fit_params'][params_names[param1_idx]]['value']
            std1 = solution['fit_params'][params_names[param1_idx]]['std']
            ax.axvline(x=mean1,linestyle='--',color='red')
            if param2_idx is not None:
                mean2 = solution['fit_params'][params_names[param2_idx]]['value']
                std2 = solution['fit_params'][params_names[param2_idx]]['std']
                ax.axhline(y=mean2,linestyle='--',color='red')
            if param2_idx is not None:
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

    if suppressAxes:
        pl.yticks([])
        pl.xticks([])

    return ax

def plot_2Ddistribution(ax, fit_out, params_names, param1_idx, param2_idx, plot_contour=False, suppressAxes=False, color='Blues',log_cmap=False):

    #old code

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
    
    #customising the colormap
    new_cmap = pl.cm.get_cmap(color)
    new_cmap.set_under('w')
    if log_cmap is True:
        cmap_norm = mplot_colors.LogNorm()
    else:
        cmap_norm = mplot_colors.Normalize()
    
    zd = grid_density_boxsum(np.min(allx), np.min(ally), np.max(allx), np.max(ally), 256, 256, zip(allx, ally))
    ax.imshow(zd, origin='lower', cmap=new_cmap,vmin=1e-10, norm=cmap_norm,
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

