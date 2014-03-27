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
import pickle
import os
from library.library_preselector import *


class preselector(object):
    def __init__(self, params):
        #
        self.params = params




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
            print '1 done'
        if generateSpectra == True:
            generate_spectra_lib(self.params,self.params.pre_cross_path,self.params.pre_speclib_path,
                                 MIXING=self.params.pre_mixing_ratios)
            print '2 done'
        if generatePCA == True:
            generate_PCA_library(self.params.pre_speclib_path+'*',self.params.pre_pca_path)
            print '3 done'


    def generate_mask(self):
        ble = True

