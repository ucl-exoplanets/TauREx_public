###############################
#Library containing plotting routines that are otherwise too cumbersome for the main code
#
###############################

import  os, sys,string, glob, itertools
import pylab as pl
import numpy as np
from numpy import array
import scipy.ndimage as ndimage
from matplotlib.ticker import FuncFormatter,ScalarFormatter
import logging

try:  
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

scale_pow = 0

def plot_multinest_results(nestdir, parameters,  save2pdf=False, out_path=None, plot_contour=True):

        # new multinest plotting routine. Inherited from plot_chains.py

        nest_raw = pymultinest.Analyzer(n_params=len(parameters),
                                        outputfiles_basename=os.path.join(nestdir, '1-'))

        nested = {}
        # n_params = nest_raw.n_params
        #only one thread for nested sampling hence nested[0]
        nested[0]= {}
        nested[0]['state'] = nest_raw.get_stats()
        nest_data = np.loadtxt(os.path.join(nestdir, '1-.txt'))#nest_raw.get_equal_weighted_posterior()
        n_params = len(nest_data[0,:])-2
        params_nest = []
        for i in range(n_params):
            nested[0][str(i)] ={}
            nested[0][str(i)]['data'] = nest_data[:,i+2]
            nested[0][str(i)]['stats'] ={}
            for j in range(len(nested[0]['state']['modes'])): #number of independent modes detected
                nested[0][str(i)]['stats'][j] ={}
                nested[0][str(i)]['stats'][j]['sigma'] = nested[0]['state']['modes'][j]['sigma'][i]
                nested[0][str(i)]['stats'][j]['mean'] = nested[0]['state']['modes'][j]['mean'][i]
                if j == 1:
                    # only take the first two modes if many are detected
                    break
            params_nest.append(str(i))

        fig_nest = plot_posteriors(nested,n_params,params_nest,parameters,plot_contour=plot_contour,alpha=0.01)

        if save2pdf:
                filename1 = os.path.join(out_path, 'nested_posteriors.pdf')
                filename2 = os.path.join(out_path, 'nested_posteriors.jpg')
                fig_nest.savefig(filename1)
                fig_nest.savefig(filename2)
                logging.info('Plot saved in %s and %s' % (filename1, filename2))

def _plot_multinest_results(DATA, parameters,  save2pdf=False, out_path=None):

    # old multinest plotting routine

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


def plot_mcmc_results(mcmcdir, parameters,  save2pdf=False, out_path=None, plot_contour=True):
    #routine plotting MCMC chains. @todo needs to be checekd for bugs

    n_params = len(parameters) #number of parameters
    
    mcmc = {}
    for thread in glob.glob(os.path.join(mcmcdir, '*')):
        print thread
        print thread[-1]
        print

        id = int(thread[-1])
        mcmc[id] = {}
                 
        with open(os.path.join(mcmcdir, 'thread_0', 'state.txt'),'r') as statefile:
            mcmc[id]['state'] = eval(statefile.read())
                 
        chainlist = glob.glob(os.path.join(thread, 'Chain_0', '*'))
        for trace in chainlist:
            traceid = string.rsplit(trace, '/')[-1][:-4]
            mcmc[id][traceid] ={}
            mcmc[id][traceid]['data'] = np.loadtxt(trace)
            mcmc[id][traceid]['stats'] = {}
            mcmc[id][traceid]['stats'][0] = {}
            mcmc[id][traceid]['stats'][0]['sigma'] = np.std(mcmc[id][traceid]['data'])
            mcmc[id][traceid]['stats'][0]['mean'] = np.mean(mcmc[id][traceid]['data'])
#                 print id, traceid, 'mean: ', np.mean(mcmc[id][traceid]['data'])
#                 print id, traceid, 'std: ', np.std(mcmc[id][traceid]['data'])

                
        
    #     print mcmc.keys()
    params_mcmc = mcmc[0]['state']['stochastics'].keys()
    params_mcmc[1:] = np.sort(params_mcmc[1:])
    
    fig_mcmc = plot_posteriors(mcmc,n_params,params_mcmc,parameters,plot_contour=plot_contour,alpha=0.05)

    if save2pdf:
            filename = os.path.join(out_path, 'mcmc_posteriors.pdf')
            fig_mcmc.savefig(filename)
            logging.info('Plot saved in %s' % filename)



def _plot_mcmc_results(DATA, parameters=False, save2pdf=False, out_path=None):
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
            pl.xlabel(parameters[i])
            pl.ylabel(parameters[j])
            for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_rotation(-30)

    if save2pdf:
            filename = os.path.join(out_path, 'mcmc_posteriors.pdf')
            fig.savefig(filename)
            logging.info('Plot saved in %s' % filename)


def plot_posteriors(data,n_params,params,display_params,plot_contour,alpha=0.05):

    fig = pl.figure(figsize=(5*n_params, 5*n_params))
#     pl.title('Nested posteriors')
    seq = 0
    globalxlims = []
    for i in range(n_params):
        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        seq +=1
#         ax.annotate('i ='+str(i)+' n ='+str(n_params * i + i + 1)+' s ='+str(seq),xy=(0.7,0.1),xycoords='axes fraction',
#                         horizontalalignment='center', verticalalignment='center')
        ax = plot_1Dposterior(ax, data, params[i],confidence=plot_contour)
    #     nested_plot.plot_marginal(i, with_ellipses=True, with_points=True, use_log_values=False, grid_points=50)
        if i == 0:
            pl.ylabel("Prob. density",fontsize=20)
        pl.xlabel(display_params[i],x=0.6, ha='right', fontsize=20)
        globalxlims= ax.get_xlim()
#         print globalxlims
       
        
        #scaling axis labels @todo labels can be improved but works now at least
        pl.rc('font', size=20)
        SciFormatter = ScalarFormatter(useMathText=True,useOffset=True,useLocale=True)
#         SciFormatter.set_scientific(True)
        SciFormatter.set_powerlimits((-1, 4))
        ax.get_xaxis().set_major_formatter(SciFormatter)
        ax.get_xaxis().get_offset_text().set_x(0.9)
#         ax.set_xlabel('{0} ({1})'.format(ax.get_xlabel(), ax.get_xaxis().get_offset_text().get_text()))
        
#         print ax.xaxis.get_offset_text().get_text()
        
#         print SciFormatter.offset
        
    #     ax.ticklabel_format(style='sci')
    #     print np.int(np.round(np.log10(mcmc[0][params[i]][0])))   
    #     print nested[0][params[i]][0]
#        print np.log10(data[0][params[i]]['data'][0])

#        if scale_pow < 0.0:
#            ax.get_xaxis().set_major_formatter(FuncFormatter(exp_formatter_fun))
#            ax.set_xlabel(display_params[i]  + ' (x $10^{{{0:d}}})$'.format(scale_pow))

#         DEPRECIATED CODE
#         scale_tmp = np.round(np.log10(data[0][params[i]]['data'][0]))
#         if np.isfinite(scale_tmp):
#             scale_pow = np.int(scale_tmp) 
#         else:
#             scale_pow = 0.0
#     #     print scale_pow
#         
#         if scale_pow < 0.0:
#             ax.get_xaxis().set_major_formatter(FuncFormatter(exp_formatter_fun))
#             ax.set_xlabel(display_params[i]  + ' (x $10^{{{0:d}}})$'.format(scale_pow))

#             ax.set_xlabel(params[i]  + ' (x $10^{{{0:d}}})$'.format(scale_pow))


        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation(+40)
            tick.label.set_fontsize(20)

        if i == 0:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_rotation(+40)
                tick.label.set_fontsize(20)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_visible(False)


        #removing overlapping tick labels
        xticks = ax.xaxis.get_major_ticks()
    #     xticks[0].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)
    #     xticks[-2].label1.set_visible(False)
    #     if i != 0:
    #         yticks = ax.yaxis.get_major_ticks()
    #         yticks.label.set_visible(False)
    #     yticks[0].label1.set_visible(False)
    #     yticks[-1].label1.set_visible(False)
    #     yticks[-2].label1.set_visible(False)


        for j in range(i):
            ax2 = pl.subplot(n_params, n_params, n_params * j + i + 1)

            seq += 1
#             ax2.annotate('i ='+str(i)+' j ='+str(j)+' n ='+str(n_params * j + i + 1)+' s ='+str(seq),xy=(0.7,0.1),xycoords='axes fraction',
#                         horizontalalignment='center', verticalalignment='center')
            ax2.ticklabel_format(style='sci')
            plot_2Ddistribution(ax2,data, 0,params[i],params[j],suppressAxes=True,confidence=plot_contour,alpha=alpha)
    #         nested_plot.plot_conditional(i, j,with_ellipses=False, with_points=True, grid_points=50)
            pl.subplots_adjust(hspace=0.0,wspace=0.0)
            ax2.set_xlim(globalxlims)

    #             plot_2Ddistribution(0,params[i],params[j],suppressAxes=True)
    #         if j ==0:
    #             pl.yticks()
    #         plot(DATA.trace(names[i])[:],DATA.trace(names[j])[:],'.',color='k')
    #         pl.xlabel(params[i])
    #         pl.ylabel(params[j])
            for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_rotation(+30)
                tick.label.set_fontsize(20)

    #         else:
    #             print i, j
    #             ax2.set_xlim(globalxlims[j])
    return fig




def plot_2Ddistribution(axis,data,thread,varname0,varname,confidence=False,suppressAxes=False,alpha=0.05):

#     std1 = data[0][varname0]['sigma']
#     mean1 = data[0][varname0]['mean']
#     std2 = data[0][varname]['sigma']
#     mean2 = data[0][varname]['mean']

    #plot 2D distributions
    axis.plot(data[thread][varname0]['data'],data[0][varname]['data'],'k.',alpha=alpha)



    #plot 1,2,3 sigma contours
    if confidence:

        for i in range(len(data[0][varname]['stats'])):
           mean1 = data[0][varname0]['stats'][i]['mean']
           mean2 = data[0][varname]['stats'][i]['mean']
           axis.axvline(x=mean1,linestyle='--',color='red')
           axis.axhline(y=mean2,linestyle='--',color='red')

        postcov = np.cov(data[0][varname0]['data'],data[0][varname]['data'])
        [eigval,eigvec] = np.linalg.eig(postcov)
        postangle = np.arctan2(eigvec[1,:],eigvec[0,:])
#         corrcoeff = np.corrcoef(mcmc[0][varname0], mcmc[0][varname])

        for i in range(len(data[0][varname]['stats'])):
            mean1 = data[0][varname0]['stats'][i]['mean']
            mean2 = data[0][varname]['stats'][i]['mean']
            std1  = data[0][varname0]['stats'][i]['sigma']
            std2  = data[0][varname]['stats'][i]['sigma']

            eX1,eY1 = ellipse(std1,std2,postangle[0],mean1,mean2)
            eX2,eY2 = ellipse(2.0*std1,2.0*std2,postangle[0],mean1,mean2)
            eX3,eY3 = ellipse(3.0*std1,3.0*std2,postangle[0],mean1,mean2)

            axis.plot(eX1,eY1,linewidth=1.0,color='red',linestyle='-')
            axis.plot(eX2,eY2,linewidth=1.0,color='orange',linestyle='-')
            axis.plot(eX3,eY3,linewidth=1.0,color='yellow',linestyle='-')

            pl.ylim([np.min(data[0][varname]['data']),np.max(data[0][varname]['data'])])
            pl.xlim([np.min(data[0][varname0]['data']),np.max(data[0][varname0]['data'])])

#         axis.annotate(corrcoeff[0,1],xy=(0.8,0.1),xycoords='axes fraction',
#                     horizontalalignment='center', verticalalignment='center')

    if suppressAxes:
        pl.yticks([])
        pl.xticks([])


def plot_1Dposterior(axis,data,varname,confidence):
#     std1 = data[0][varname]['sigma']
#     mean1 = data[0][varname]['mean']
    for i in data.keys():
        axis.hist(data[i][varname]['data'],bins=50,normed=True,alpha=0.8, linewidth=2.0,histtype='step') # @todo custom bins?
#         pl.title(varname)
    if confidence:
        for i in range(len(data[0][varname]['stats'])):
    #         print data[0][varname]['stats'][i]
            mean1 = data[0][varname]['stats'][i]['mean']
            std1 = data[0][varname]['stats'][i]['sigma']
            axis.axvline(x=mean1,linestyle='--',color='red')
    return axis

#@todo the plot label bug must be somewhere around here
def exp_formatter_fun(x, p):

    global scale_pow

    if x < 1.0:
        return "%.2f" % (x * (10 ** np.abs(scale_pow)))
#        try:
#            scale = np.int(np.round(np.log10(x)))
#        except OverflowError:
#            scale = x
#        print x, (10 ** np.abs(scale)), (x * (10 ** np.abs(scale)))

#        return "%.2f" % (x * (10 ** np.abs(scale)))
#             return "%.2f" % (x / np.exp(scale_pow))
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



def iterate_TP_profile(fit_params, fit_params_std, fit_idx, TP_function):
    '''
    function iterating through all lower and upper bounds of parameters
    to determine which combination gives the lowest/highest attainable 
    TP profile. Returns mean TP profile with errorbars on each pressure level
    ''' 
    
    TP_params     = fit_params[fit_idx:]
    TP_params_std = fit_params_std[fit_idx:]
    
    
    Tmean,P,X = TP_function(fit_params)
    
    bounds = [] #list of lower and upper parameter bounds 
    for i in xrange(len(TP_params)):
        bounds.append((TP_params[i]-TP_params_std,TP_params[i]+TP_params_std))
        
        
    iterlist = list(itertools.product(*bounds))
    iter_num = np.shape(iterlist)[0] #number of possible combinations
    
    T_iter   = np.zeros((len(Tmean),iter_num))
    T_minmax = np.zeros((len(Tmean),2))
    
    for i in range(iter_num):
        T,P, X = TP_function(iterlist[i])
        T_iter[:,i] = T
    
    T_minmax[:,0] = np.min(T_iter)
    T_minmax[:,1] = np.max(T_iter)
    
    return Tmean, T_minmax, P
    
    
    
    
    
    




