################################################
#class preselector
#
# This class has multiple functions:
# - it correlates data to library of pre-computed single molecule spectra
#   to estimate number of molecules and rough abundances as priors to
#   later analysis
# - it generates lookup libraries for correlation analysis
# - it handles all required pre-processing steps required by main code
#
# Input: -parameter object
#
#
#
# Output: - Priors for main code (to be defined)
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Feb 2014
#
################################################

#loading libraries
import numpy as np
import pylab as pl
import pickle,gzip,os
from library.library_preselector import *


class preselector(object):
    def __init__(self, params,data):
        #
        self.params = params
        self.data   = data

        #reading in spectrum data to be fitted
        self.spectrum = data.spectrum
        self.nwave = len(self.spectrum[:,0])
        self.wavegrid = self.spectrum[:,0]




    def run_preprocess(self,convertLinelist=None,generateSpectra=None,generatePCA=None):
        #function running pre-processing steps
        #convertLinelist will convert ExoMol cross sections from wavenumbers to microns
        #generateSpectra will generate library of transmission spectra from linelists
        #generatePCA will run PCA on spectral library and compile final dictionary

        #getting booleans from parameter file if not set manually
        if convertLinelist == None:
            convertLinelist = self.params.pre_conver2microns
        if generateSpectra == None:
            generateSpectra = self.params.pre_gen_speclib
        if generatePCA == None:
            generatePCA = self.params.pre_gen_pca

        #doing the pre-processing
        if convertLinelist == True:
            convert2microns(self.params.pre_cross_path+'*')
            # print '1 done'
        if generateSpectra == True:
            generate_spectra_lib(self.params,self.params.pre_cross_path,self.params.pre_speclib_path,
                                 MIXING=self.params.pre_mixing_ratios)
            # print '2 done'
        if generatePCA == True:
            generate_PCA_library(self.params.pre_speclib_path+'*',self.params.pre_pca_path)
            # print '3 done'


    def load_library(self,PATH=None):
        #function loading spec_pcalib.pkl.zip into memory
        if PATH == None:
            PATH = self.params.pre_pca_path

        try:
            with gzip.open(PATH+'spec_pcalib.pkl.zip',mode='rb') as filehandle:
                self.PCALIB = pickle.load(filehandle)
        except IOError:
            print 'WARNING: cannot find spec_pcalib.pkl.zip in: ',PATH
            redo = raw_input('Try to generate library from scratch? (y/n) [n]: ')
            if redo == 'N' or redo == 'n' or redo == 'no' or redo == '':
                exit()
            elif redo == 'Y' or redo == 'y' or redo == 'yes':
                print 'generating library...'
                self.run_preprocess(convertLinelist=True,generateSpectra=True,generatePCA=True)
                print 'loading library...'
                with gzip.open(PATH+'spec_pcalib.pkl.zip',mode='rb') as filehandle:
                    self.PCALIB = pickle.load(filehandle)


    def interpolate2data(self):
        #interpolates the normalised PCs of PCALIB to the observed data grid
        molkeys = self.PCALIB.keys()
        for molecule in molkeys:
            normpc     = self.PCALIB[molecule]['PCA']['norm']
            pcwavegrid = self.PCALIB[molecule]['wavegrid']
            interpc    = np.zeros((self.nwave,np.shape(normpc)[1]))
            for i in range(np.shape(normpc)[1]):
                interpc[:,i] = np.interp(self.wavegrid,pcwavegrid,normpc[:,i])
            self.PCALIB[molecule]['PCA']['norm_interp'] = interpc



    def generate_mask(self):
        #generates mask for the correlation processs
        #the masks are generated from the first normalised PC

        molkeys = self.PCALIB.keys()
        for molecule in molkeys:
            pc1 = self.PCALIB[molecule]['PCA']['norm'][:,0]
            pc2 = self.PCALIB[molecule]['PCA']['norm'][:,1]

            mask = pc1 > 0.6

            print mask

            # print pc1
            # exit()
            pl.figure(1)
            pl.plot(pc1)
            pl.plot(pc2)

            pl.figure(2)
            pl.plot(pc1[mask])
            pl.show()



