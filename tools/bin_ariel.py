#small program to bin your taurex output to the ariel wavelength grid and plot the emission spectrum and TP-profiles 

import numpy as np 
import pylab as pl
import cPickle as pkl
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
    def __init__(self,FILEPATH,FILENAME,NOISEFILE,SNR=None):
       #initialising bin_spectrum class 
        
        #setting constants 
        self.RJUP = 6.9911e7
        self.RSOL = 6.955e8
        self.MJUP = 1.898e27
        self.G = 6.67384e-11
        self.KBOLTZ = 1.380648813e-23

        
        #reading in instance file 
        print 'reading in spectrum: {}'.format(FILENAME)
        self.FILEIO = os.path.join(FILEPATH,FILENAME)
        self.read_instance_file(self.FILEIO)
        
        self.fn_type = self.data['params']['gen_type']
        self.cmap = cm.get_cmap('Paired')
        
        print 'reading in noise file: {}'.format(NOISEFILE)
        self.NOISEIO = os.path.join(FILEPATH,NOISEFILE)
        self.read_noise_model(self.NOISEIO)
        
        #generating ariel wavelength grid 
        print 'generating ARIEL wavelength grid...'
        self.build_ariel_lambda_grid()
        
        #binning down internal high resolution spectrum to ariel grid 
        print 'binning spectrum to ARIEL resolution...'
        self.highres_spec = self.data['data']['spectrum']
        self.lowres_spec = self.bin_spectrum(self.highres_spec, self.lambda_grid[:,0],self.lambda_grid[:,1])
        
        print 'interpolating noise to wavelength grid...'
        self.interpolate_noise()
                
        if SNR is not None:
            print 'scaling noise to SNR: {}'.format(SNR)
            if self.fn_type == 'emission':
                self.scale_noise_emission(np.float(SNR))
            elif self.fn_type == 'transmission':
                self.scale_noise_transmission(np.float(SNR))
        
        #putting everything together into array: 
        #wavelength, flux, noise, wavelength_bin_size
        self.final = self.create_output()
        
    
    
    #internal class functions 
    def read_instance_file(self,FILENAME):
        with open(FILENAME,'r') as f:
            self.data = pkl.load(f)
    
    def read_noise_model(self,NOISE_FILE):
        #Reading in ENZO noise csv file
        csvlist =[]
        with open(NOISE_FILE) as f:
            csvreader = csv.reader(f)
            for row in csvreader:
                csvlist.append(row)
        
        csvlist = csvlist[7:]
        csvlist2 = []
        for x in csvlist:
            if x[0] is '':
                x[0] = '-1'
            csvlist2.append(list(map(float,x)))
        csvarr = np.asarray(csvlist2)
            
        self.raw_wave = csvarr[:,1]
        self.noise_phot = csvarr[:,2]
        self.noise_all = csvarr[:,3]
    
    
    def build_ariel_lambda_grid(self):
        #function generating the ariel wavelenght grid 
        # 1.25-1.95 micron R = 10
        # 1.95-3.9 micron R = 100
        # 3.9-7.8 micron R = 30
        
        R1 = 10
        R2 = 100
        R3 = 30
        
        [val1_l,idx1_l] = self.find_nearest(self.raw_wave,1.25)
        [val1_u,idx1_u] = self.find_nearest(self.raw_wave,1.95)
        
        [val2_l,idx2_l] = self.find_nearest(self.raw_wave,1.95)
        [val2_u,idx2_u] = self.find_nearest(self.raw_wave,3.9)
        
        [val3_l,idx3_l] = self.find_nearest(self.raw_wave,3.9)
        [val3_u,idx3_u] = self.find_nearest(self.raw_wave,7.8)
        
        #first spectral band 
        lam_1,dlam_1 = self.fill_lam_dlam(val1_l,val1_u,R1)
        #second spectral band
        lam_2,dlam_2 = self.fill_lam_dlam(val2_l,val2_u,R2)
        #third spectral band
        lam_3,dlam_3 = self.fill_lam_dlam(val3_l,val3_u,R3)
        #photommetric bands
        lam_phot = [0.525,0.9,1.125]
        dlam_phot = [0.05,0.2,0.15]
        
        #getting sizes of bands 
        S0 = len(lam_phot)
        S1 = len(lam_1)
        S2 = len(lam_2)
        S3 = len(lam_3)
        
        #building output array and filling it 
        out_band = np.zeros(((S0+S1+S2+S3),2))
        out_band[:S0,0] = lam_phot
        out_band[:S0,1] = dlam_phot
        
        out_band[S0:S1+S0,0] = lam_1
        out_band[S0:S1+S0,1] = dlam_1
        
        out_band[S0+S1:S0+S1+S2,0] = lam_2
        out_band[S0+S1:S0+S1+S2,1] = dlam_2
        
        out_band[S0+S1+S2:,0] = lam_3
        out_band[S0+S1+S2:,1] = dlam_3
        
        #sorting grid 
        specgrid_idx = np.argsort(out_band[:,0])
        self.lambda_grid = out_band[specgrid_idx,:]
    
    def interpolate_noise(self):
        #takes Enzo's noise file for the planet and interpolates this onto the ARIEL wavelength grid 
        self.noise_int_all = np.interp(self.lambda_grid[:,0],self.raw_wave,self.noise_all)
        self.noise = self.noise_int_all
    
    def get_atm_scale_height(self):
        #calculates the atmospheric scaleheight
        self.planet_radius = self.data['params']['planet_radius'] 
        self.planet_mass   = self.data['params']['planet_mass'] 

        self.g = (self.G * self.planet_mass) / (self.planet_radius**2) # surface gravity (0th layer)
        self.H = (self.KBOLTZ*np.mean(self.data['data']['temperature_profile'][:,1]))/(self.data['data']['mu_profile'][0,1]*self.g) # scaleheight at the surface (0th layer)
    
    
    def get_mean_noise_flux(self):
        #function calculating the mean nnoise and mean flux accross a band
        
        #finding the indices for the spectroscopic bands
        [val1_l,idx1_l] = self.find_nearest(self.raw_wave,1.25)
        [val1_u,idx1_u] = self.find_nearest(self.raw_wave,1.95)
        
        [val2_l,idx2_l] = self.find_nearest(self.lambda_grid[:,0],1.95)
        [val2_u,idx2_u] = self.find_nearest(self.lambda_grid[:,0],3.9)
        
        [val3_l,idx3_l] = self.find_nearest(self.lambda_grid[:,0],3.9)
        [val3_u,idx3_u] = self.find_nearest(self.lambda_grid[:,0],7.8)
        
        #calculating average noise over each band2 (1.95-3.9 microns) and band2 (3.9 - 7.8microns)
        noise_band2_aver = np.mean(self.noise[idx2_l:idx2_u])
        noise_band3_aver = np.mean(self.noise[idx3_l:idx3_u])
        #calculating average signal over each band2 (1.95-3.9 microns) and band2 (3.9 - 7.8microns)
        signal_band2_aver = np.mean(self.lowres_spec[idx2_l:idx2_u,1])
        signal_band3_aver = np.mean(self.lowres_spec[idx3_l:idx3_u,1])
        
        return noise_band2_aver, noise_band3_aver, signal_band2_aver,signal_band3_aver
        
    
    
    def scale_noise_transmission(self,SNR_final=20):
        #function scaling the noise model to a certain signal-to-noise per band. 
        #the signal is defined as 5 x scale_height of planet 
        
        #getting planetary scaleheight 
        self.get_atm_scale_height()
        
        #getting mean noise and flux accross bands
        noise_band2_aver, noise_band3_aver, signal_band2_aver,signal_band3_aver = self.get_mean_noise_flux()
        
        Rs = self.data['params']['star_radius']
        RpH = self.planet_radius + 5*self.H
        #calculating the atmosperhic fraction
        ATMfrac = (RpH**2 - self.planet_radius**2) / Rs**2

        SNR2 = ATMfrac/noise_band2_aver
        SNR3 = ATMfrac/noise_band3_aver
        
        #determining where the SNR is higher, band2 or band3 and calculating a scaling factor to scale the 
        #lowest SNR band to the required SNR_final
        if SNR2 >= SNR3:
            scale_fact = SNR_final/SNR3
        else:
            scale_fact = SNR_final/SNR2
 
        #scaling entire noise array to SNR_final
        self.noise /= scale_fact
    
    def scale_noise_emission(self,SNR_final=20):
        #function scaling the noise model to a certain signal-to-noise per band. 
        #the SNR in each band is defined as 
        #SNR_band = mean(flux_in_band)/error
        
        #getting mean noise and flux accross bands
        noise_band2_aver, noise_band3_aver, signal_band2_aver,signal_band3_aver = self.get_mean_noise_flux()

        #calculating average SNR
        SNR2 = signal_band2_aver/noise_band2_aver
        SNR3 = signal_band3_aver/noise_band3_aver
        
        #determining where the SNR is higher, band2 or band3 and calculating a scaling factor to scale the 
        #lowest SNR band to the required SNR_final
        if SNR2 >= SNR3:
            scale_fact = SNR_final/SNR3
        else:
            scale_fact = SNR_final/SNR2
    
#         print 'scale ',scale_fact
#         noise_band2_aver /= scale_fact
#         noise_band3_aver /= scale_fact
        
        #scaling entire noise array to SNR_final
        self.noise /= scale_fact
        


    def create_output(self):
        #creating new array: wavelength, flux, noise, wavelength_bin_size
        out = np.zeros((len(self.lambda_grid[:,0]),4))
        out[:,0] = self.lowres_spec[:,0]
        out[:,1] = self.lowres_spec[:,1]
        out[:,2] = self.noise
        out[:,3] = self.lambda_grid[:,1]
        
        return out
        
            
    def plot_spectrum(self,plot_errors=False):
        #plotting spectrum 
        pl.figure()
        if plot_errors: 
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c='k',label='high-res',alpha=0.1)
            (_, caps, _) = pl.errorbar(self.final[:,0],self.final[:,1],yerr=self.final[:,2],c='b',linewidth=2.5, elinewidth=3.0,label='ARIEL-res')
#             pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k.',linewidth=2.5,markersize=10.0)   
            for cap in caps:
                cap.set_markeredgewidth(2)
        else:
            pl.plot(self.highres_spec[:,0],self.highres_spec[:,1],c=[0.8,0.8,0.8],label='high-res')
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k',linewidth=2.5,label='ARIEL-res')
            pl.plot(self.lowres_spec[:,0],self.lowres_spec[:,1],'k.',linewidth=2.5,markersize=10.0)
        pl.legend(loc=0)
        pl.xlim([0.5,8.0])
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
        
    def save_spectrum(self,OUTPATH,OUTFILE='SPECTRUM_ARIEL.dat'):
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


    def bin_spectrum(self,spectrum, wlgrid,dlam_grid):
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
    parser.add_argument('--noise',
                      dest='noise_file',
                      default='Ariel noise on transit - MCR v2 - HD189733b.csv',
                       help='Enzo\'s noise file (.csv)'
                      ) 
    parser.add_argument('--dir',
                      dest='instance_dir',
                      default='',
                      help='input file directory'
                      )
    parser.add_argument('--snr',
                      dest='snr',
                      default=None,
                      help='minimum SNR in bands 2 & 3'
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
                      help='saves final spectrum to ascii')
    parser.add_argument('--outfile',
                        dest='out_file',
                        default='SPECTRUM_ARIEL.dat',
                        help='output filename')
    parser.add_argument('--outdir',
                        dest='out_dir',
                        default='',
                        help='output directory')
    
    #parsing command line options
    options = parser.parse_args()
    
    
    #initialising binning object 
    bin_ob = bin_spectrum(FILEPATH=options.instance_dir,
                          FILENAME=options.instance_file,
                          NOISEFILE=options.noise_file,
                          SNR=options.snr)
    
    
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
    
    
    
    