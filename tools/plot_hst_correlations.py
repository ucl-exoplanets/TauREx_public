#! /usr/bin/python

import numpy as np
import pylab as pl
import pickle as pkl
import string as str
import itertools as iter
import mpl_toolkits.mplot3d.axes3d as axes3d
from mpl_toolkits.mplot3d import proj3d
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import os 
import glob
import argparse
from __builtin__ import False




class hst_correlations(object):
    def __init__(self,FILEPATH='./',SNR=None):
       #initialising bin_spectrum class 
        
        #setting constants 
        self.RJUP = 6.9911e7
        self.RSOL = 6.955e8
        self.MJUP = 1.898e27
        self.G = 6.67384e-11
        self.KBOLTZ = 1.380648813e-23
        
        #setting variables 
        self.filepath = FILEPATH
        plot_params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20,10),
         'axes.labelsize': '15',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'15',
         'ytick.labelsize':'15'}
        pl.rcParams.update(plot_params)
        self.labels = options.labels
        
        colors = [(1, 0, 0), (0, 1, 0)]  # R -> G 
        n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
        cmap_name = 'red_green'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)
        self.cmap=cm
        
        #loading data 
        self.load_data()
        
        #calculating derived parameters 
        self.get_derived_params()
      
#         print self.data[self.planet_list[0]]['derived'].keys()
#         exit()
      
    
    def get_derived_params(self):
        #calculates derived parameters
        
        for i in range(self.N_instance):
            paramlist = self.data[self.planet_list[i]]['fit_params_names']
            paramlist.append('mu_derived')
            
            tracedata = self.data[self.planet_list[i]]['solutions'][0]['tracedata']
            tracedata = np.column_stack((tracedata, self.data[self.planet_list[i]]['solutions'][0]['fit_params']['mu_derived']['trace']))
            weights   = self.data[self.planet_list[i]]['solutions'][0]['weights']
    
            self.data[self.planet_list[i]]['derived'] = {}
            comb = iter.combinations(paramlist,2)
            for per in comb:
                name = '{0}-{1}'.format(per[0],per[1])
                self.data[self.planet_list[i]]['derived'][name] = {} 
                #pearson correlation coefficient 
                corr_tmp = self.weight_corr(tracedata[:,paramlist.index(per[0])],tracedata[:,paramlist.index(per[1])],weights[:])
                self.data[self.planet_list[i]]['derived'][name]['pearson'] = corr_tmp 
                #linear regression 
                fit_tmp = self.fit_linear_weighted(tracedata[:,paramlist.index(per[0])],tracedata[:,paramlist.index(per[1])], weights[:])
                self.data[self.planet_list[i]]['derived'][name]['regression'] = fit_tmp[1]
                
    def fit_linear_weighted(self,data1,data2,weights):
        #fits linear regression line to weighted data
        fit = np.polynomial.polynomial.polyfit(data1, data2, deg=1, w=weights)
        return fit
    
    def cov(self,x, y, w):
    #weighted covariance
        return np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(x, weights=w))) / np.sum(w)
    
    def weight_corr(self,x, y, w):
     #weighted pearson corrletaion coefficient 
        return self.cov(x, y, w) / np.sqrt(self.cov(x, x, w) * self.cov(y, y, w))
    
    def get_fitted_list(self):
        #prints out fitted parameters
        print 'Fitted parameters:'
        print self.data[self.planet_list[0]]['solutions'][0]['fit_params'].keys()
        
    def get_static_list(self):
        #prints out static parameters 
        print 'Static parameters:'
        print self.data[self.planet_list[0]]['params'].keys()
        
    def get_derived_list(self):
        #prints out static parameters 
        print 'Derived parameters:'
        print self.data[self.planet_list[0]]['derived'].keys()
    
    
    
    def plot_scatter(self,options,colour='k',save_plot=False):
        #plots either 2D or 3D scatter plots given number of parameters 
        
        PARAMS = options.static
        FIT_PARAMS = options.fitted
        DERIVED_TYPE = options.derived_type
        DERIVED = options.derived
        
        if len(FIT_PARAMS) > 0 and FIT_PARAMS[0] == 'logE':
            data,parlist = self.compile_data(PARAMS=['planet_mass'], FIT_PARAMS=[],DERIVED_TYPE=None,DERIVED_PARAMS=[])
            fig = self.plot_logE(data)
            if save_plot:
                pl.savefig('{0}.pdf'.format('logEvidence'))
        else:
            data,parlist = self.compile_data(PARAMS, FIT_PARAMS,DERIVED_TYPE,DERIVED)
            if len(parlist) == 2:
                fig = self.plot_scatter_1D(data,parlist,colour)
                if save_plot:
                    pl.savefig('{0}.pdf'.format(parlist[0]))
            if len(parlist) == 4: 
                fig = self.plot_scatter_2D(data,parlist,colour)
                if save_plot:
                    pl.savefig('scatter-{0}-{1}.pdf'.format(parlist[0],parlist[2]))
            elif len(parlist) == 6:
                fig = self.plot_scatter_3D(data,parlist,colour)
                if save_plot:
                    pl.savefig('scatter-{0}-{1}-{2}.pdf'.format(parlist[0],parlist[2],parlist[4]))
            else: 
                print 'Wrong number of paramters, has to be either 2 or 3!'
        
        
    def plot_logE(self,data):
        #plots goodness of fit 
        fig = pl.figure()
        xax = np.asarray([i for i in range(self.N_instance)])
        pl.bar(xax,data[:,-1],color='g')
        pl.xticks(xax,self.planet_list,rotation=60)
        pl.xlabel('Planet name')
        pl.ylabel('log(Evidence)')
        
    def plot_scatter_1D(self,data,parlist,colour):
        #plotting 2D scatter with errorbars
        norm_logE = (data[:,-1]-np.min(data[:,-1]))/np.max(data[:,-1]-np.min(data[:,-1]))

        fig = pl.figure()
        xax = np.asarray([i for i in range(self.N_instance)])
        pl.errorbar(xax,data[:,0],yerr=data[:,1],marker='o',linestyle='',linewidth=0.0)
        for i in range(self.N_instance):
            pl.errorbar(xax[i],data[i,0],yerr=data[i,1],marker='o',linestyle='',color=self.cmap(norm_logE[i]))
            cax = pl.plot(xax[i],data[i,0],color=self.cmap(norm_logE[i]),marker='o',markersize=10)
#         cax = pl.plot(xax,data[:,0],color=self.cmap(norm_logE[i]),marker='o',markersize=10)


        pl.xticks(xax,self.planet_list,rotation=60)
        
        pl.xlabel('Planet name')
        pl.ylabel(parlist[0])
        pl.title('log(Evidence): {0} - {1}'.format(np.min(data[:,-1]),np.max(data[:,-1])))
#         pl.colorbar(cax,[np.min(data[:,-1]),np.max(data[:,-1])])
    
        return fig
    
    def plot_scatter_2D(self,data,parlist,colour):
        #plotting 2D scatter with errorbars
        norm_logE = (data[:,-1]-np.min(data[:,-1]))/np.max(data[:,-1]-np.min(data[:,-1]))

        fig = pl.figure()
        
        for i in range(self.N_instance):
            pl.errorbar(data[i,0],data[i,2],xerr=data[i,1],yerr=data[i,3],marker='o',linestyle='',color=self.cmap(norm_logE[i]))
            cax = pl.plot(data[i,0],data[i,2],color=self.cmap(norm_logE[i]),marker='o',markersize=10)
        pl.xlabel(parlist[0])
        pl.ylabel(parlist[2])
        pl.title('log(Evidence): {0} - {1}'.format(np.min(data[:,-1]),np.max(data[:,-1])))
        
        if self.labels:
            for label, x, y in zip(self.planet_list, data[:, 0], data[:, 2]):
                pl.annotate(
                    label, 
                    xy = (x, y), xytext = (-20, 20),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                    arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        
#         pl.colorbar(cax,[np.min(data[:,-1]),np.max(data[:,-1])])
        return fig
        
    def plot_scatter_3D(self,data,parlist,colour='k'):
        #plotting 3D scatter with errorbars
        norm_logE = (data[:,-1]-np.min(data[:,-1]))/np.max(data[:,-1]-np.min(data[:,-1]))
        
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel(parlist[0])
        ax.set_ylabel(parlist[2])
        ax.set_zlabel(parlist[4])
        ax.set_title('log(Evidence): {0} - {1}'.format(np.min(data[:,-1]),np.max(data[:,-1])))
        
        fx = []; fy = []; fz = [];
        xerror = []; yerror = []; zerror = [];
        
        for i in range(self.N_instance):
            fx.append(data[i,0])
            xerror.append(data[i,1])
            fy.append(data[i,2])
            yerror.append(data[i,3])
            fz.append(data[i,4])
            zerror.append(data[i,5])
        
        #plot errorbars
        for i in np.arange(0, len(fx)):
            ax.plot([fx[i]+xerror[i], fx[i]-xerror[i]], [fy[i], fy[i]], [fz[i], fz[i]], marker="_",color=self.cmap(norm_logE[i]))
            ax.plot([fx[i], fx[i]], [fy[i]+yerror[i], fy[i]-yerror[i]], [fz[i], fz[i]], marker="_",color=self.cmap(norm_logE[i]))
            ax.plot([fx[i], fx[i]], [fy[i], fy[i]], [fz[i]+zerror[i], fz[i]-zerror[i]], marker="_",color=self.cmap(norm_logE[i]))
            ax.plot([fx[i], fx[i]], [fy[i], fy[i]], [fz[i], fz[i]],color=self.cmap(norm_logE[i]),marker='o',markersize=10,linestyle="None")
            
        if self.labels:
            for i in range(len(fx)):
                text=self.planet_list[i]    
                x2, y2, _ = proj3d.proj_transform(fx[i],fy[i],fz[i], ax.get_proj())    
                label = pl.annotate(text,
                xycoords='data', 
                xy = (x2, y2), xytext = (60, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        return fig
        
        
    def compile_data(self,PARAMS,FIT_PARAMS,DERIVED_TYPE=None,DERIVED_PARAMS=None):
        #routine that compiles data subsets from input keywords 
        #array format:
        #parameter1, sig_parameter1,..., fitted1,sig_fitted1,..., logE
        N_params    = len(PARAMS) * 2 #value - error pairs
        N_fit       = len(FIT_PARAMS) * 2 #value - error pairs
        N_totalpar  = N_fit + N_params+1
        if DERIVED_TYPE is not None:
            N_derived = len(DERIVED_PARAMS) * 2 
            N_totalpar += N_derived
        else: 
            DERIVED_PARAMS = []

        
        out = np.zeros((self.N_instance,N_totalpar))
        parlist_out = []

        
        idx1 = 0
        for planet in self.planet_list:
            idx2 = 0
            for param in PARAMS:
                out[idx1,idx2] = self.data[planet]['params'][param]
                out[idx1,idx2+1] = 0.0
                if idx1 == 0:
                    parlist_out.append(param)
                    parlist_out.append(param+'_sig')
                idx2 += 2
            for fit in FIT_PARAMS:
                out[idx1,idx2] = self.data[planet]['solutions'][0]['fit_params'][fit]['nest_mean']
                out[idx1,idx2+1] = self.data[planet]['solutions'][0]['fit_params'][fit]['nest_sigma']
                if idx1 == 0:
                    parlist_out.append(fit)
                    parlist_out.append(fit+'_sig')
                idx2 += 2
            
            for derived in DERIVED_PARAMS:
                out[idx1,idx2] = self.data[planet]['derived'][derived][DERIVED_TYPE]
                out[idx1,idx2+1] = 0.0
                if idx1 == 0:
                    parlist_out.append(derived)
                    parlist_out.append(derived+'_sig')
                idx2 += 2
                
                
            out[idx1,-1] = self.data[planet]['global_logE'][0]
            idx1 += 1
        
        return out,parlist_out
        
    def load_data(self):
        #main routine reading in all nest_out.pickle files and loading them into the data dictionary 
        self.instance_list = glob.glob(os.path.join(self.filepath,'*'))
        self.N_instance = 0
        self.planet_list = []
        
        #setting up main data dictionary. THis contains the raw nest_out.pickle dictionaries
        #from which later the relevant arrays will be compiled. 
        self.data = {} 
        for dir in self.instance_list: 
            planet_tmp = str.split(dir,'/')[-1]
            planet = str.split(planet_tmp,'_')[0]
            self.data[planet] = {}
            data = self._read_instance_file(os.path.join(dir,'1'))
            if data is 'False':
                data = self._read_instance_file(os.path.join(dir,'0'))
            if data is not 'False':
                self.data[planet] = data
                self.planet_list.append(planet)
                self.N_instance += 1
            
            
        
    def _read_instance_file(self,FILEPATH,FILENAME='nest_out.pickle'):
        #subroutine reading in nest_out.pickle files 
        try: 
            with open(os.path.join(FILEPATH,FILENAME),'r') as f:
                return pkl.load(f)
        except IOError:
            print 'CANNOT LOAD: ',os.path.join(FILEPATH,FILENAME)
            return 'False'
        
        
            


if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()


    parser.add_argument('--dir',
                      dest='instance_dir',
                      default='',
                      help='input file directory'
                      )
    parser.add_argument('--static',
                      dest='static',
                      default=[],
                      nargs='+',
                      help='list of static parameters to be plotted'
                      )
    parser.add_argument('--fitted',
                      dest='fitted',
                      default=[],
                      nargs='+',
                      help='list of fitted parameters to be plotted'
                      )
    parser.add_argument('--derived',
                      dest='derived',
                      default=[],
                      nargs='+',
                      help='list of derived parameters to be plotted'
                      )
    parser.add_argument('--d_type',
                        dest='derived_type',
                        default='regression',
                        help='type of derived parameter: regression, pearson. Default: regression'
                        )
    parser.add_argument('--list_static',
                        dest='list_static',
                        action='store_true',
                        default=False,
                        help='list all static parameters'
                        )
    parser.add_argument('--list_fitted',
                        dest='list_fitted',
                        action='store_true',
                        default=False,
                        help='list all fitted parameters'
                        )
    parser.add_argument('--list_derived',
                        dest='list_derived',
                        action='store_true',
                        default=False,
                        help='list all derived parameters'
                        )
    
    parser.add_argument('--save_plot',
                      action='store_true',
                      dest='save_plot',
                      default=False,
                      help='saves plot to pdf'
                      )
    parser.add_argument('--labels',
                      action='store_true',
                      dest='labels',
                      default=False,
                      help='put labels on figure'
                      )
#     parser.add_argument('--save',
#                       action='store_true',
#                       dest='save_file',
#                       default=False,
#                       help='saves final spectrum to ascii')
#     parser.add_argument('--outfile',
#                         dest='out_file',
#                         default='SPECTRUM_ARIEL.dat',
#                         help='output filename')
#     parser.add_argument('--outdir',
#                         dest='out_dir',
#                         default='',
#                         help='output directory')
    
    #parsing command line options
    options = parser.parse_args()
    
#     print options.static
#     exit()
    
    
    #initialising binning object 
    corr_ob = hst_correlations(FILEPATH=options.instance_dir)      
    
    #returning parameter names 
    if options.list_fitted:
        corr_ob.get_fitted_list()
    if options.list_static:
        corr_ob.get_static_list()
    if options.list_derived:
        corr_ob.get_derived_list()
        
    #plotting stuff
    if len(options.fitted) > 0 or len(options.static) > 0 or len(options.derived) > 0:
        print 'Close window to terminate program...'
        fig = corr_ob.plot_scatter(options)
        pl.show()
    
          