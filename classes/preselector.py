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


class preselector(object):
    def __init__(self, params):
        #question on how autonomous this class should be... should most pre-processing be done
        #on start up or on explicit demand? Needs some form of history to avoid unnecessary duplication.
        self.params = params

        self.LIB = self.load_spectrum_library
        self.save_spectrum_library()


    def load_spectrum_library(self):
        #checks if spectrum library exists if it does it is loaded. if not it is created
        if os.path.isfile(self.params.in_pre_file):
            with open(self.params.in_pre_file, 'rb') as LIBfilehandle:
                LIB = pickle.load(LIBfilehandle)
            print 'loaded from file'
        else:
            LIB = {}
            LIB.history = {}
            LIB.type = {}
            LIB.type.trans = {}
            LIB.type.trans.mol = {}
            LIB.type.emis = {}
            LIB.type.emis.mol = {}
            print 'created library'
        return LIB


    def save_spectrum_library(self):
        #saving library to file
        with open(self.params.in_pre_file, 'wb') as LIBfilehandle:
            pickle.dump(self.LIB, LIBfilehandle)


    def check_history(self):
        #checking if preprocessing steps have already been exectured in past
        ble = 1

    def generate_priors(self):
        #collects priors from correlation analysis and passes them to main code in standardised format
        ble = 1

    def correlate(self):
        ble = 1

    def create_masks(self):
        #search for 'important' features in molecule line list by thresholding (?!)
        #create wavelength grid of important features per molecule to be used for correlation analysis
        ble = 1

    def generate_library_grid(self):
        #genrate temperature, mixing ratio grid for available molecules
        #generate indexing map of library (dictionaries or arrays?)
        ble = 1

    def generate_library(self):
        #load preexisting
        #load grid and input ABS files
        #save generated to pickled object
        #save generated to ascii files

        ble = 1
