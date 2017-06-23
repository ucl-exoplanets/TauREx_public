#script plotting some things from the create_spectrum.py instance 


import numpy as np 
import pylab as pl
import cPickle as pkl
import os 
import argparse
import csv
from matplotlib import cm
from matplotlib import rc
import matplotlib as mpl
from __builtin__ import False
from matplotlib.ticker import ScalarFormatter


#some global matplotlib vars
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['text.antialiased'] = True
#rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 11})



class plot_spectrum(object):
    def __init__(self,FILEPATH,FILENAME,BINSPEC=False,BINRES=300,SPEC_LOG=False):
       #initialising bin_spectrum class 
        
        self.bin_spec = BINSPEC
        self.spec_log = SPEC_LOG
        #setting constants 
        self.RJUP = 6.9911e7
        self.RSOL = 6.955e8
        self.MJUP = 1.898e27
        self.G = 6.67384e-11
        self.KBOLTZ = 1.380648813e-23
        self.cmap = cm.get_cmap('Paired')
        self.cmap2 = cm.get_cmap('gray')
        self.BINRES = BINRES
        
        
        
        #reading in instance file 
        print 'reading in spectrum: {}'.format(FILENAME)
        self.FILEIO = os.path.join(FILEPATH,FILENAME)
        self.read_instance_file(self.FILEIO)
        self.fn_type = self.data['params']['gen_type']
        self.highres_spec = self.data['data']['spectrum'][::-1]
        
        #setting xaxis bounds 
        self.xmin = 0.4
        self.xmax = 25.0
        if np.min(self.highres_spec[:,0]) > self.xmin:
            self.xmin = np.min(self.highres_spec[:,0])
        if np.max(self.highres_spec[:,0]) < self.xmax:
            self.xmax = np.max(self.highres_spec[:,0])
        
        
        #binning down internal high resolution spectrum to ariel grid 
        if self.bin_spec:
            print 'binning spectrum to resolution: {0}'.format(BINRES)
            self.get_specgrid(BINRES, self.highres_spec[0,0], self.highres_spec[-1,0])
            self.lowres_spec = self.bin_spectrum(self.highres_spec, self.lambda_grid)

        print 'creating output grid...'
        self.final = self.create_output()    
    
    #internal class functions 
    def read_instance_file(self,FILENAME):
        with open(FILENAME,'r') as f:
            self.data = pkl.load(f)
            
    def get_atm_scale_height(self):
        #calculates the atmospheric scaleheight
        self.planet_radius = self.data['params']['planet_radius'] 
        self.planet_mass   = self.data['params']['planet_mass'] 

        self.g = (self.G * self.planet_mass) / (self.planet_radius**2) # surface gravity (0th layer)
        self.H = (self.KBOLTZ*np.mean(self.data['data']['temperature_profile'][:,1]))/(self.data['data']['mu_profile'][0,1]*self.g) # scaleheight at the surface (0th layer)


    def create_output(self):
        #creating new array: wavelength, flux, noise, wavelength_bin_size
        if self.bin_spec:
            out = np.zeros((len(self.lambda_grid),2))
            out[:,0] = self.lowres_spec[:,0]
            out[:,1] = self.lowres_spec[:,1]
        else:
            out = np.zeros((len(self.highres_spec[:,0]),2))
            out[:,0] = self.highres_spec[:,0]
            out[:,1] = self.highres_spec[:,1]
        return out
        
    def get_specgrid(self,R=5000, lambda_min=0.1, lambda_max=20.0):
        #generating wavelength grid with uniform binning in log(lambda)
        #lambda min and max are in microns, R is the spectral resolution
        #R = lambda/delta_lamdba
        specgrid = []
        delta_lambda =[]
        specgrid.append(lambda_min)
        run = True
        i=1
        while run:
            dlam= specgrid[i-1]/np.float(R)
            specgrid.append(specgrid[i-1]+dlam)
            delta_lambda.append(dlam)
    
            if specgrid[i] >= lambda_max:
                run=False
            i+=1

        self.lambda_grid = np.asarray(specgrid)
        self.dlambda_grid = np.asarray(delta_lambda)

            
    def plot_spectrum(self,plot_errors=False):
        #plotting spectrum 
        pl.figure()
        if self.bin_spec:
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c=[0.8,0.8,0.8],label='high-res')
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k',linewidth=2.5,label='low-res')
#             pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k.',linewidth=2.5,markersize=10.0)
        else:
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],'k',label='high-res')
        pl.legend(loc=0)
      
        pl.xlim([self.xmin,self.xmax])
        pl.title('Spectrum')
        if self.spec_log:
            pl.xscale('log')
        if self.fn_type == 'emission':
            pl.ylabel(r'$F_p/F_\ast$')
        else:
            pl.ylabel(r'$(R_p/R_s)^2$')
        pl.xlabel(r'Wavelength ($\mu$m)')
        
    def plot_contributions(self):
        #plotting contributions 
        N = len(self.data["data"]["opacity_contrib"])
        
        fig = pl.figure(figsize=(7,3.5))
        ax = fig.add_subplot(111)

#         plot_spectrum = np.zeros((len(self.highres_spec[:,0]), 2))
#         plot_spectrum[:,0] = self.highres_spec[:,0]
#         plot_spectrum[:,1] = self.highres_spec[:,1]
#         plot_spectrum = plot_spectrum[plot_spectrum[:,0].argsort(axis=0)] # sort in wavelength
#         if self.data['params']['in_opacity_method'][:4] == 'xsec':
#             self.get_specgrid(self.BINRES, self.highres_spec[0,0], self.highres_spec[-1,0])
#             self.lowres_spec = self.bin_spectrum(self.highres_spec, self.lambda_grid)
# #             plot_spectrum = binspectrum(plot_spectrum, self.BINRES)
#         #plt.plot(plot_spectrum[:,0], plot_spectrum[:,1], alpha=0.7, color='black', label='Fitted model')
        if self.bin_spec:
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k',linewidth=2.5,label='low-res')
        else:
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c=[0.8,0.8,0.8],label='high-res')
#         plt.errorbar(obs[:,0], obs[:,1], obs[:,2], lw=1, color='black', alpha=0.5, ls='none', zorder=99, label='Observed')

        for idx, val in enumerate(self.data["data"]["opacity_contrib"]):
            plot_spectrum = self.data["data"]["opacity_contrib"][val]
            plot_spectrum = plot_spectrum[plot_spectrum[:,0].argsort(axis=0)] # sort in wavelength
            if self.bin_spec:
                plot_spectrum = self.bin_spectrum(plot_spectrum, self.lambda_grid)
#                 plot_spectrum = binspectrum(plot_spectrum, 300)
            pl.plot(plot_spectrum[:,0], plot_spectrum[:,1], color=self.cmap(float(idx)/N), label=val)

#         pl.xlim(np.min(obs[:,0])-0.05*np.min(obs[:,0]), np.max(obs[:,0])+0.05*np.max(obs[:,0]))
        pl.xlabel('Wavelength ($\mu$m)')
        pl.ylabel('$(R_p/R_*)^2$')
        pl.tick_params(axis='x', which='minor')
        # set log scale only if interval is greater than 5 micron
        if self.spec_log:
            pl.xscale('log')
            pl.tick_params(axis='x', which='minor')
#             ax.xaxis.set_minor_formatter(FormatStrFormatter("%i"))
#             ax.xaxis.set_major_formatter(FormatStrFormatter("%i"))
        pl.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, prop={'size':11}, frameon=False)
        pl.xlim([self.xmin,self.xmax])
    
        
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
        
    def plot_transmittance(self):
        
        #getting wavelength grid 
        wave = self.data['data']['spectrum'][:,0]
        if wave[0] > wave[-1]:
            wave = wave[::-1]
        
        #getting pressure grid 
        pressure_alt = self.data['data']['altitude_profile']
        pressure = pressure_alt[:,0][::-1]
#         print(pressure)
        
        #getting transmittance grid
        transmit = self.data['data']['transmittance']
        
        #dumping to file 
        np.savetxt('transmittance_wave.dat',wave)
        np.savetxt('transmittance_pressure.dat',pressure)
        np.savetxt('transmittance.dat',transmit)
        
        
        fig = pl.figure(figsize=(10,5))
        ax = fig.add_subplot(111)
        mesh = ax.pcolormesh(wave, pressure/1e2,transmit, cmap=cm.coolwarm)
        
#         pl.xlim(0.3, 50)
#         pl.ylim(100, 1e-7)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Wavelength ($\mu$m)')
        ax.set_ylabel('Pressure (mbar)')
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.invert_yaxis()
#         pl.gca().annotate('Hot Jupiter spectral transmittance', xy=(.05, .94),
#                 xycoords='axes fraction', color='white',
#                 horizontalalignment='left', verticalalignment='top',
#                 fontsize=12)
        
        ax.set_xlim([self.xmin,self.xmax])
        cbar = pl.colorbar(mesh,ax=ax)
        cbar.set_label('Transmittance')
#         pl.savefig('Transmittance.pdf', format='pdf')
        
#         pl.gca().xaxis.set_ticks([0.5, 1, 5, 10, 20 , 30, 40 , 50])
#         pl.gca().xaxis.set_ticklabels(['0.5', '1', '5', '10', '20', '30', '40', '50'])
#         pl.gca().invert_yaxis()
#         pl.yscale('log')
        
        
    def save_spectrum(self,OUTPATH,OUTFILE='SPECTRUM_ARIEL.dat'):
        #saving the spectrum to some output file 
        outname = os.path.join(OUTPATH,OUTFILE)
        np.savetxt(outname,self.final)     

    def find_nearest(self,arr, value):
        # find nearest value in array
        arr = np.array(arr)
        idx = (abs(arr-value)).argmin()
        return [arr[idx], idx]
    
    
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
            for i in range(len(wavegrid)-1):
                bingrid.append(wavegrid[i]-binwidths[i]/2.)
                bingrid.append(wavegrid[i]+binwidths[i]/2.)
    
            # build bin grid index array (an index for each model datapoint
            bingrid_idx = np.empty(len(specgrid))
            bingrid_idx[:] = np.NaN
            for i in range(len(specgrid)):
                for j in range(len(wavegrid)-1):
                    if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                        bingrid_idx[i] = j+1
                        break
                    
#         print bingrid_idx
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
    parser.add_argument('--bin',
                      dest='bin_spec',
                      action='store_true',
                      default=False,
                       help='bin spectrum to differnt resolution'
                      ) 
    parser.add_argument('--res',
                        dest='resolution',
                        default=100,
                        help='Spectral resolution of binning (constant in delta_lambda)'
                        )
    parser.add_argument('--spec_log',
                        dest='spec_log',
                        default=False,
                        action='store_true',
                        help='plots spectrum in log(wavelength)')
    parser.add_argument('--dir',
                      dest='instance_dir',
                      default='',
                      help='input file directory'
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
    plot_ob = plot_spectrum(FILEPATH=options.instance_dir,
                          FILENAME=options.instance_file,
                          BINSPEC=options.bin_spec,
                          BINRES=options.resolution,
                          SPEC_LOG=options.spec_log)
    
    
    #saving final spectrum to file
    if options.save_file:
        print 'saving file to: {}'.format(os.path.join(options.out_dir,options.out_file))
        plot_ob.save_spectrum(options.out_dir, options.out_file)
        
    #plotting spectrum and tp-profile
  
    print 'plotting stuff...'
    plot_ob.plot_spectrum(plot_errors=True)
    plot_ob.plot_tp_profile()
    plot_ob.plot_xprofiles()
    try:
        plot_ob.plot_contributions()
    except KeyError:
        pass
    try: 
        plot_ob.plot_transmittance()
    except KeyError:
        pass
    
    pl.show()