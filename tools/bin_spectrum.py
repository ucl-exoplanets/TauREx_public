#! /usr/bin/pyton

#script spitting out model comparisions between nemesis and taurex for emission

import numpy as np 
import pylab as pl
import cPickle as pkl 
import subprocess as sub
import glob
import os

def get_specgrid( R=5000, lambda_min=0.1, lambda_max=20.0):
    #generating wavelength grid with uniform binning in log(lambda)
    #lambda min and max are in microns, R is the spectral resolution
    #R = lambda/delta_lamdba
    specgrid = []
    delta_lambda =[]
    specgrid.append(lambda_min)
    run = True
    i=1
    while run:
        dlam= specgrid[i-1]/R
        specgrid.append(specgrid[i-1]+dlam)
        delta_lambda.append(dlam)

        if specgrid[i] >= lambda_max:
            run=False
        i+=1
    return np.asarray(specgrid),np.asarray(delta_lambda)

def get_specbingrid(wavegrid, specgrid, binwidths=None):
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
        # this bingrid is actually useless, as it doesn't allow for gaps in the data
        bingrid = []
        for i in range(len(wavegrid)):
            bingrid.append(wavegrid[i]-binwidths[i]/2.)
            bingrid.append(wavegrid[i]+binwidths[i]/2.)

        # build bin grid index array (an index for each model datapoint
        bingrid_idx = np.empty(len(specgrid))
        print 'BLEEE ',len(specgrid)
        bingrid_idx[:] = np.NaN
        for i in range(len(specgrid)):
            for j in range(len(wavegrid)):
                if specgrid[i] >= (wavegrid[j]-binwidths[j]/2.) and specgrid[i] < (wavegrid[j]+binwidths[j]/2.):
                    bingrid_idx[i] = j+1
                    break
    return bingrid, bingrid_idx


def bin_spectrum(spectrum, wlgrid, noise=0):
    sp_wlgrid = spectrum[:,0]
    bin_up =  (sp_wlgrid[-1]-sp_wlgrid[-2])/2.
    bin_low = (sp_wlgrid[1]-sp_wlgrid[0])/2.
    lambdamin = sp_wlgrid[0] - bin_low
    lambdamax = sp_wlgrid[-1] + bin_up
    sp_bingrid, sp_bingrididx = get_specbingrid(wlgrid, sp_wlgrid)
    sp_nbingrid = len(wlgrid)
    model_binned = [spectrum[:,1][sp_bingrididx == i].mean() for i in xrange(1, sp_nbingrid+1)]
    if noise > 0:
        model_binned += np.random.normal(0, noise, len(model_binned))
    return model_binned



def binspectrum(spectrum_in, resolution):
    wavegrid, dlamb_grid = get_specgrid(R=resolution,lambda_min=np.min(spectrum_in[:,0]),lambda_max=np.max(spectrum_in[:,0]))
    spec_bin_grid, spec_bin_grid_idx = get_specbingrid(wavegrid, spectrum_in[:,0])
    spectrum_binned = [spectrum_in[:,1][spec_bin_grid_idx == i].mean() for i in xrange(1,len(spec_bin_grid))]
    return transpose(vstack((spec_bin_grid[:-1], spectrum_binned)))


BINRES = 100

DIR = '/Users/ingowaldmann/Desktop/lighttemp/'

data_list = glob.glob(DIR+'*.dat')

# data_list = [DIR+'Modern-Earth-Lightning_NO.dat', DIR+'Modern-Earth-Lightning_NO2.dat', DIR+'Modern-Earth-Lightning_CO.dat', DIR+'Modern-Earth-Lightning_CO2.dat']
data_list = [DIR+'Modern-Earth-Lightning_O2.dat']


print data_list

# exit()

for specfile in data_list:
    
    spectrum = np.loadtxt(specfile)
    
    specgrid,delta_lambda = get_specgrid(R=BINRES,lambda_min=np.min(spectrum[:,0]),lambda_max=np.max(spectrum[:,0]))
    model = bin_spectrum(spectrum, specgrid, noise=0.0)
    
    out = np.zeros((len(specgrid), 2))
    out[:,0] = specgrid
    out[:,1] = model
    out = out[~np.isnan(out).any(1)]
    
    np.savetxt(specfile[:-4]+'_R'+str(BINRES)+'.dat',out)
    
#     pl.figure(1)
#     pl.plot(model)
#     pl.show()
