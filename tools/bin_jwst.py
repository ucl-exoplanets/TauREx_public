#small program to bin your taurex output to the jwst wavelength grid and plot the spectrum and TP-profiles 

import numpy as np 
import pylab as pl
import cPickle as pkl
import pyfits as pf 
import os 
import argparse
import csv
from matplotlib import cm
from matplotlib import rc
import matplotlib as mpl

#some global matplotlib vars
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['text.antialiased'] = True
#rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 11})


class bin_spectrum(object):
    def __init__(self,FILEPATH,FILENAME,NOISEFILE):
       #initialising bin_spectrum class 
        
        #setting constants 
        self.RJUP = 6.9911e7
        self.RSOL = 6.955e8
        self.MJUP = 1.898e27
        self.G = 6.67384e-11
        self.KBOLTZ = 1.380648813e-23
        self.cmap = cm.get_cmap('Paired')
        
        #reading in instance file 
        print 'reading in spectrum: {}'.format(FILENAME)
        self.FILEIO = os.path.join(FILEPATH,FILENAME)
        self.read_instance_file(self.FILEIO)
        self.fn_type = self.data['params']['gen_type']
        
        print 'reading in JWST file: {}'.format(NOISEFILE)
        self.NOISEIO = os.path.join(FILEPATH,NOISEFILE)
        self.read_jwst_fits(self.NOISEIO)
        
        #binning down internal high resolution spectrum to ariel grid 
        print 'binning spectrum to JWST resolution...'
        self.highres_spec = self.data['data']['spectrum']
        self.lowres_spec = self.bin_spectrum(self.highres_spec, self.lambda_grid)
        
        print 'creating output grid...'
        self.final = self.create_output()    
    
    #internal class functions 
    def read_instance_file(self,FILENAME):
        with open(FILENAME,'r') as f:
            self.data = pkl.load(f)
            
    def read_jwst_fits(self,FNAME):
        #loading PO-Lagage provided .fits files 
        fdata = pf.open(FNAME)
        data = fdata[1].data[0] #reads main data table
        gidx = [data[1] != 0] #excludes all zero flux regions
        self.lambda_grid  = data[0][gidx] #wavelength array
        self.noise_jwst   = data[1][gidx] #noise array 
        
    def get_atm_scale_height(self):
        #calculates the atmospheric scaleheight
        self.planet_radius = self.data['params']['planet_radius'] 
        self.planet_mass   = self.data['params']['planet_mass'] 

        self.g = (self.G * self.planet_mass) / (self.planet_radius**2) # surface gravity (0th layer)
        self.H = (self.KBOLTZ*np.mean(self.data['data']['temperature_profile'][:,1]))/(self.data['data']['mu_profile'][0,1]*self.g) # scaleheight at the surface (0th layer)


    def create_output(self):
        #creating new array: wavelength, flux, noise, wavelength_bin_size
        out = np.zeros((len(self.lambda_grid),3))
        out[:,0] = self.lowres_spec[:,0]
        out[:,1] = self.lowres_spec[:,1]
        out[:,2] = self.noise_jwst
        
        return out
        
            
    def plot_spectrum(self,plot_errors=False):
        #plotting spectrum 
        pl.figure()
        if plot_errors: 
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c='k',label='high-res',alpha=0.1)
            (_, caps, _) = pl.errorbar(self.final[:,0],self.final[:,1],yerr=self.final[:,2],c='b',linewidth=2.5, elinewidth=3.0,label='JWST-res')
#             pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k.',linewidth=2.5,markersize=10.0)   
            for cap in caps:
                cap.set_markeredgewidth(2)
        else:
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c=[0.8,0.8,0.8],label='high-res')
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k',linewidth=2.5,label='JWST-res')
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k.',linewidth=2.5,markersize=10.0)
        pl.legend(loc=0)
        pl.xlim([0.5,15.0])
        pl.title('Spectrum')
        if self.fn_type == 'emission':
            pl.ylabel(r'$F_p/F_\ast$')
        else:
            pl.ylabel(r'$(R_p/R_s)^2$')
        pl.xlabel(r'Wavelength ($\mu$m)')
        
    def plot_tp_profile(self):
        #plotting temperature-pressure profile 
        tp = self.data['data']['temperature_profile']
        pl.figure()
        pl.plot(tp[:,1],tp[:,0]/100000.0,linewidth=2.5)
        pl.xlim([np.min(tp[:,1])-200.0,np.max(tp[:,1])+200.0])
        pl.gca().invert_yaxis()
        pl.yscale('log')
        pl.xlabel('Temperature (K)')
        pl.ylabel('Pressure (bar)')
        
    def plot_xprofiles(self):
        #plotting abundance profiles 
        fig = pl.figure()
        ax = fig.add_subplot(111)
        N = len(self.data['params']['atm_active_gases']+self.data['params']['atm_inactive_gases'])
        cols_mol = {}
        for mol_idx, mol_val in enumerate(self.data['params']['atm_active_gases']):
            cols_mol[mol_val] = self.cmap(float(mol_idx)/N)
            prof = self.data['data']['active_mixratio_profile'][mol_idx]
            if mol_val in cols_mol:
                pl.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], label=mol_val, ls='solid',linewidth=2.0)
            else:
                pl.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(len(self.data['params']['atm_active_gases'])+len(self.data['params']['atm_inactive_gases'])+mol_idx)/N), label=mol_val, ls='solid',linewidth=2.0)

        for mol_idx, mol_val in enumerate(self.data['params']['atm_inactive_gases']):
            prof = self.data['data']['inactive_mixratio_profile'][mol_idx]
            if mol_val in cols_mol:
                pl.plot(prof[:,1], prof[:,0]/1e5, color=cols_mol[mol_val], ls='dashed',linewidth=2.0)
            else:
                pl.plot(prof[:,1], prof[:,0]/1e5, color=self.cmap(float(len(self.data['params']['atm_active_gases'])+len(self.data['params']['atm_inactive_gases'])+mol_idx)/N), label=mol_val, ls='dashed',linewidth=2.0)

        pl.gca().invert_yaxis()
        pl.yscale('log')
        pl.xscale('log')
        pl.xlim(1e-12, 3)
        pl.xlabel('Mixing ratio')
        pl.ylabel('Pressure (bar)')
        pl.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'size':11}, frameon=False)
        pl.title('Mixing ratios', fontsize=14)
        
    def save_spectrum(self,OUTPATH,OUTFILE='SPECTRUM_JWST.dat'):
        #saving the spectrum to some output file 
        outname = os.path.join(OUTPATH,OUTFILE)
        np.savetxt(outname,self.final)
        
    
    def find_nearest(self,arr, value):
        # find nearest value in array
        arr = np.array(arr)
        idx = (abs(arr-value)).argmin()
        return [arr[idx], idx]
    
    def fill_lam_dlam(self,val_low,val_up,res):
        #version constant in resolution
        lam = []
        dlam = []
        lam.append(val_low)
        dlam.append(val_low/res) 
    
        run = True
        while run:
            dtmp = lam[-1]/res
            if (lam[-1] + dtmp) > val_up:
                break
    
            lam.append(lam[-1]+dtmp)
            dlam.append(dtmp)
        return lam,dlam
    
    def fill_lam_dlam_static(self,val_low,val_up,res):
        #version constant in dlam
        lam = []
        dlam = []
        lam.append(val_low)
        dtmp = ((val_low+val_up)/2.0)/res 
        dlam.append(dtmp)
        
        run = True 
        while run:
            if (lam[-1] + dtmp) > val_up:
                break
            lam.append(lam[-1]+dtmp)
            dlam.append(dtmp)
            
        return lam, dlam
    
    def get_specbingrid(self,wavegrid, specgrid, binwidths=None):
        #function calculating the bin boundaries for the data
        #this is used to bin the internal spectrum to the data in fitting module
        if not isinstance(binwidths, (np.ndarray, np.generic)):
            bingrid =[]
            bingrid.append(wavegrid[0]- (wavegrid[1]-wavegrid[0])/2.0) #first bin edge
            for i in range(len(wavegrid)-1):
                bingrid.append(wavegrid[i]+(wavegrid[i+1]-wavegrid[i])/2.0)
            bingrid.append((wavegrid[-1]-wavegrid[-2])/2.0 + wavegrid[-1]) #last bin edge
            bingrid_idx = np.digitize(specgrid,bingrid) #getting the specgrid indexes for bins
        else:
            bingrid = []
            for i in range(len(wavegrid)):
                bingrid.append(wavegrid[i]-binwidths[i]/2.)
                bingrid.append(wavegrid[i]+binwidths[i]/2.)
    
            # build bin grid index array (an index for each model datapoint
            bingrid_idx = np.empty(len(specgrid))
    #         print 'BLEEE ',len(specgrid)
            bingrid_idx[:] = np.NaN
            for i in range(len(specgrid)):
                for j in range(len(wavegrid)):
                    if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                        bingrid_idx[i] = j+1
                        break
        return bingrid, bingrid_idx


    def bin_spectrum(self,spectrum, wlgrid,dlam_grid=None):
        sp_wlgrid = spectrum[:,0]
        sp_bingrid, sp_bingrididx = self.get_specbingrid(wlgrid, sp_wlgrid,dlam_grid)
        sp_nbingrid = len(wlgrid)
        model_binned = [spectrum[:,1][sp_bingrididx == i].mean() for i in xrange(1, sp_nbingrid+1)]
        
        out = np.zeros((len(wlgrid),2))
        out[:,0] = wlgrid
        out[:,1] = model_binned
        return out




if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--instance',
                      dest='instance_file',
                      default='SPECTRUM_INSTANCE_out.pickle',
                       help='TauREx create_spectrum instance file'
                      )
    parser.add_argument('--jwst',
                      dest='jwst_file',
                      default='JWST noise file',
                       help='PO-Lagage\'s noise file (.fits)'
                      ) 
    parser.add_argument('--dir',
                      dest='instance_dir',
                      default='',
                      help='input file directory'
                      )
    parser.add_argument('--plot',
                      action='store_true',
                      dest='plot',
                      default=False,
                      help='turns on plotting'
                      )
    parser.add_argument('--save',
                      action='store_true',
                      dest='save_file',
                      default=False,
                      help='saves final spectrum to ascii'
                      )
    parser.add_argument('--outfile',
                        dest='out_file',
                        default='SPECTRUM_ARIEL.dat',
                        help='output filename'
                        )
    parser.add_argument('--outdir',
                        dest='out_dir',
                        default='',
                        help='output directory'
                        )
    
    #parsing command line options
    options = parser.parse_args()
    
    
    #initialising binning object 
    bin_ob = bin_spectrum(FILEPATH=options.instance_dir,
                          FILENAME=options.instance_file,
                          NOISEFILE=options.jwst_file)
    
    
    #saving final spectrum to file
    if options.save_file:
        print 'saving file to: {}'.format(os.path.join(options.out_dir,options.out_file))
        bin_ob.save_spectrum(options.out_dir, options.out_file)
        
    #plotting spectrum and tp-profile
    if options.plot:
        print 'plotting stuff...'
        bin_ob.plot_spectrum(plot_errors=True)
        bin_ob.plot_tp_profile()
        bin_ob.plot_xprofiles()
        pl.show()