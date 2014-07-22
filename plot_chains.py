#! /usr/bin/python 

import glob, string
import numpy as np
import pylab as pl
from numpy import array
import scipy.ndimage as ndimage
from matplotlib.ticker import FuncFormatter


#function time
# Define a function to make the ellipses
def ellipse(ra,rb,ang,x0,y0,Nb=100):
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=np.linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y



def plot_2Ddistribution(axis,data,thread,varname0,varname,confidence=False,suppressAxes=False):
    
#     std1 = data[0][varname0]['sigma']
#     mean1 = data[0][varname0]['mean']
#     std2 = data[0][varname]['sigma']
#     mean2 = data[0][varname]['mean']
    
    #plot 2D distributions
    axis.plot(data[thread][varname0]['data'],data[0][varname]['data'],'k.',alpha=0.05)
    
    for i in range(len(data[0][varname]['stats'])):
        mean1 = data[0][varname0]['stats'][i]['mean']
        mean2 = data[0][varname]['stats'][i]['mean']
        axis.axvline(x=mean1,linestyle='--',color='red')
        axis.axhline(y=mean2,linestyle='--',color='red')

    #plot 1,2,3 sigma contours
    if confidence:
       
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
        

def plot_1Dposterior(axis,data,varname):
#     std1 = data[0][varname]['sigma']
#     mean1 = data[0][varname]['mean']
    for i in data.keys():
        axis.hist(data[i][varname]['data'],bins=BINS,normed=True,alpha=0.8, linewidth=2.0,histtype='step')
#         pl.title(varname)
    for i in range(len(data[0][varname]['stats'])):
#         print data[0][varname]['stats'][i]
        mean1 = data[0][varname]['stats'][i]['mean']
        axis.axvline(x=mean1,linestyle='--',color='red')
    return axis
    
    
def exp_formatter_fun(x, p):
        if x < 1.0:
            try:
                scale = np.int(np.round(np.log10(x))) 
            except OverflowError:
                scale = x
#             print x, (10 ** np.abs(scale)), (x * (10 ** np.abs(scale)))
            return "%.2f" % (x * (10 ** np.abs(scale)))
#             return "%.2f" % (x / np.exp(scale_pow))
        if x >= 1.0:
            return "%.2f" % (x)
        
        
        
def plot_posteriors(data,n_params,params,display_params):
    fig = pl.figure(figsize=(5*n_params, 5*n_params))
#     pl.title('Nested posteriors')
    seq = 0
    globalxlims = []
    for i in range(n_params):
        ax = pl.subplot(n_params, n_params, n_params * i + i + 1)
        seq +=1
#         ax.annotate('i ='+str(i)+' n ='+str(n_params * i + i + 1)+' s ='+str(seq),xy=(0.7,0.1),xycoords='axes fraction',
#                         horizontalalignment='center', verticalalignment='center')
        ax = plot_1Dposterior(ax, data, params[i])
    #     nested_plot.plot_marginal(i, with_ellipses=True, with_points=True, use_log_values=False, grid_points=50)
        if i == 0:
            pl.ylabel("Prob. density",fontsize=20)
        pl.xlabel(display_params[i],fontsize=20)
        globalxlims= ax.get_xlim()
       
        
        #scaling axis labels 
    #     ax.xaxis.get_major_formatter().set_powerlimits((0, 100))
    #     ax.ticklabel_format(style='sci')
    #     print np.int(np.round(np.log10(mcmc[0][params[i]][0])))   
    #     print nested[0][params[i]][0]
#        print np.log10(data[0][params[i]]['data'][0])
        scale_pow = np.int(np.round(np.log10(data[0][params[i]]['data'][0]))) 
    #     print scale_pow
        if scale_pow < 0.0:
            ax.get_xaxis().set_major_formatter(FuncFormatter(exp_formatter_fun))
            ax.set_xlabel(display_params[i]  + ' (x $10^{{{0:d}}})$'.format(scale_pow))
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
            plot_2Ddistribution(ax2,data, 0,params[i],params[j],suppressAxes=True,confidence=True)
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

#################################################################################


chaindir = raw_input('Chain directory [chains/]: ')
if chaindir == '':
    chaindir = 'chains/'
if chaindir[-1] != '/':
    chaindir = chaindir+'/'

mcmcdir = chaindir+'MCMC/'

plottmp = raw_input('Plot MCMC posteriors [y]/n: ')
if plottmp == '' or plottmp == 'y' or plottmp == 'Y':
    plot_mcmc = True
else:
    plot_mcmc = False
     
plottmp = raw_input('Plot Nested posteriors [y]/n: ')
if plottmp == '' or plottmp == 'y' or plottmp == 'Y':
    plot_nest = True
else:
    plot_nest = False

# chaindir = 'chains_10000/'
# plot_mcmc = False
# chaindir = 'chains/'
# mcmcdir  = 'chains/MCMC/'


# nested_post_sep = np.loadtxt(chaindir+'1-post_separate.dat')
# nested          = np.loadtxt(chaindir+'1-.txt')

BINS= 50
# n_params = 6

disp_params = np.loadtxt(chaindir+'parameters.dat',dtype='string')
disp_params[0] = disp_params[0]+' (K)' #add units to temperature parameter
n_params = len(disp_params)


if plot_mcmc:
#reading in MCMC data
    mcmc = {}
    for thread in glob.glob(mcmcdir+'*'):
        id = int(thread[-1])
        mcmc[id] = {}
             
        with open(mcmcdir+'thread_0/state.txt','r') as statefile:
            mcmc[id]['state'] = eval(statefile.read())
             
        chainlist = glob.glob(thread+'/Chain_0/*')
        for trace in chainlist:
            traceid = string.rsplit(trace, '/')[-1][:-4]
            mcmc[id][traceid] ={}
            mcmc[id][traceid]['data'] = np.loadtxt(trace)
            mcmc[id][traceid]['stats'] = {}
            mcmc[id][traceid]['stats'][0] = {}
            mcmc[id][traceid]['stats'][0]['sigma'] = np.std(mcmc[id][traceid]['data'])
            mcmc[id][traceid]['stats'][0]['mean'] = np.mean(mcmc[id][traceid]['data'])
            
    
#     print mcmc.keys()
    params_mcmc = mcmc[0]['state']['stochastics'].keys()
    params_mcmc[1:] = np.sort(params_mcmc[1:])
#     print params_mcmc
    
#     exit()
     #     mcmcchainlist = glob.glob(thread+'/Chain_0/*')[1]
     #     print string.rsplit(mcmcchainlist, '/')[-1][:-4]



# exit()


if plot_nest:
    import pymultinest as pymn
    #reading in Nested data
    nest_raw = pymn.Analyzer(n_params=n_params,outputfiles_basename=chaindir+'1-')
    nested = {}
    # n_params = nest_raw.n_params 
    #only one thread for nested sampling hence nested[0]
    nested[0]= {}
    nested[0]['state'] = nest_raw.get_stats()
    nest_data = np.loadtxt(chaindir+'1-post_equal_weights.dat')#nest_raw.get_equal_weighted_posterior()
    n_params = len(nest_data[0,:])-1
    params_nest = []
    for i in range(n_params):
        nested[0][str(i)] ={}
        nested[0][str(i)]['data'] = nest_data[:,i]
        nested[0][str(i)]['stats'] ={}
        for j in range(len(nested[0]['state']['modes'])): #number of independent modes detected 
            nested[0][str(i)]['stats'][j] ={}
            nested[0][str(i)]['stats'][j]['sigma'] = nested[0]['state']['modes'][j]['sigma'][i]
            nested[0][str(i)]['stats'][j]['mean'] = nested[0]['state']['modes'][j]['mean'][i]
            if j == 1:
                break #only take the first two modes if many are detected
        params_nest.append(str(i))
        
#     print n_params, len(params)


# pl.figure()
# pl.scatter(nested.get_equal_weighted_posterior()[:,1],nested.get_equal_weighted_posterior()[:,2])
# pl.show()
# 
# print nested[0]['0']['stats']#.keys()

# exit()
# print '-------------------------------'
# 
# print mcmc[0]['state']

# exit()

#################################################################################


# params = mcmc[0]['state']['stochastics'].keys()
# n_params = len(params)
# n_params = 3

if plot_nest: 
    fig_nest = plot_posteriors(nested,n_params,params_nest,disp_params)
    fig_nest.savefig('nested_posterior.pdf')
    fig_nest.savefig('nested_posterior.jpg')
if plot_mcmc: 
    fig_mcmc = plot_posteriors(mcmc,n_params,params_mcmc,disp_params)
    fig_mcmc.savefig('mcmc_posterior.pdf')
    fig_mcmc.savefig('mcmc_posterior.jpg')

# corrarray = np.zeros((n_params,len(mcmc[0][params[0]])))
# for i in range(n_params):
#     corrarray[i,:] = mcmc[0][params[i]]
# corrmatrix = np.corrcoef(corrarray)
# 
# # pl.figure()
# ax = pl.subplot(n_params, n_params,n_params*(n_params-1)+1)
# ax.imshow(corrmatrix,interpolation='none',cmap='gray_r')
# # ax.colorbar()

#################################################################################




# pl.show()
