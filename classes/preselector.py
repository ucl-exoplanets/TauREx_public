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
import gzip,os
import cPickle as pickle
from copy import deepcopy
from library.library_preselector import *
import scipy.stats.stats as st


class preselector(object):
    def __init__(self, params,data):
        #
        self.params = params
        self.data   = data

        #reading in spectrum data to be fitted
        self.spectrum = data.spectrum
        self.nwave = len(self.spectrum[:,0])
        self.wavegrid = self.spectrum[:,0]


    def run(self,runpreprocess=True):
        #script running the preselector

        if runpreprocess:
            #running pre_processing routines using parameter file specifications
            self.run_preprocess()
        #loading generated spectral library and running PCA
        self.load_library()
        #interpolating PCA vectors onto data vector
        self.interpolate2data()
        #generating mask for correlation using first principal component
        self.generate_mask()
        #computing euclidian distance and correlation coefficients between library and data
        self.correlate()
        #rank molecules according to their lowest euclidian distance to data
        #also estimate number of likely molecules in data
        self.rank_molecules()
        #calculating some planetary parameters
        self.calc_astroparams()

        # print self.mol_rank
        # print self.mol_idx



    def run_preprocess(self,convertLinelist=None,generateSpectra=None,generatePCA=None):
        #function running pre-processing steps
        #convertLinelist will convert ExoMol cross sections from wavenumbers to microns
        #generateSpectra will generate library of transmission spectra from linelists
        #generatePCA will run PCA on spectral library and compile final dictionary

        #getting booleans from parameter file if not set manually
        if convertLinelist is None:
            convertLinelist = self.params.pre_conver2microns
        if generateSpectra is None:
            generateSpectra = self.params.pre_gen_speclib
        if generatePCA is None:
            generatePCA = self.params.pre_gen_pca

        #doing the pre-processing
        if convertLinelist:
            convert2microns(self.params.pre_cross_path+'*')
            # print '1 done'
        if generateSpectra:
            generate_spectra_lib(self.params,self.params.pre_cross_path,self.params.pre_speclib_path,
                                 MIXING=self.params.pre_mixing_ratios)
            # print '2 done'
        if generatePCA:
            generate_PCA_library(self.params.pre_speclib_path+'*',self.params.pre_pca_path)
            # print '3 done'


    def load_library(self,PATH=None):
        #function loading spec_pcalib.pkl.zip into memory
        if PATH is None:
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



    def generate_mask(self,thres=0.5):
        #generates mask for the correlation processs
        #the masks are generated from the first normalised PC

        molkeys = self.PCALIB.keys()
        for molecule in molkeys:
            pc1 = self.PCALIB[molecule]['PCA']['norm'][:,0]

            #adding mask to object
            self.PCALIB[molecule]['PCA']['mask'] = pc1 >thres

            try:
                pc1_interp = self.PCALIB[molecule]['PCA']['norm_interp'][:,0]
                self.PCALIB[molecule]['PCA']['interp_mask'] = pc1_interp > thres
            except(KeyError):
                pass


            # print mask
            # print molecule
            # mask = pc1 > thres
            # pc2 = self.PCALIB[molecule]['PCA']['norm'][:,1]
            # xaxis = np.zeros((len(pc1),1))
            # for j in range(len(pc1)):
            #     xaxis[j] = j
            # # # print pc1
            # # # exit()
            # pl.figure(1)
            # pl.plot(xaxis,pc1,'b')
            # pl.plot(xaxis,pc2,'g')
            # pl.plot(xaxis[mask],pc1[mask],'r')
            # #
            # # # pl.figure(2)
            # # # pl.plot(pc1[mask])
            # pl.show()
            # # exit()

    def correlate(self):
        #applies pre-determined masks and correlates the data to the second
        #principal component of the spectral library

        #normalise data
        data = self.spectrum[:,1]
        datanorm = data - np.min(data)
        datanorm /= np.max(datanorm)


        molkeys = self.PCALIB.keys()
        for molecule in molkeys:
            #adding layer to PCALIB
            self.PCALIB[molecule]['preselect'] = {}

            #loading mask
            mask = self.PCALIB[molecule]['PCA']['interp_mask']

            datanorm_m = datanorm[mask]-np.mean(datanorm[mask])
            # datanorm_m /= np.max(datanorm[mask])
            pc2 = self.PCALIB[molecule]['PCA']['norm_interp'][mask,1] - np.mean(self.PCALIB[molecule]['PCA']['norm_interp'][mask,1])
            # pc2 /= np.max(pc2)
            pc2_inv = (-1.0*(pc2-np.mean(pc2)))+np.mean(pc2)

            eucdist = np.sum(sqrt((datanorm_m-pc2)**2))/len(datanorm[mask])
            eucdist_inv = np.sum(sqrt((datanorm_m-pc2_inv)**2))/len(datanorm[mask])

            corrcoeff = st.pearsonr(datanorm_m,pc2)


            self.PCALIB[molecule]['preselect']['pearson'] = corrcoeff
            # self.PCALIB[molecule]['preselect']['euclid_dist'] = eucdist
            if corrcoeff[0] <0.0:
            # if eucdist_inv < eucdist:
                self.PCALIB[molecule]['preselect']['euclid_dist'] = eucdist_inv
            else:
                self.PCALIB[molecule]['preselect']['euclid_dist'] = eucdist



            # print molecule,': ',corrcoeff, '... ',eucdist,'... ',eucdist_inv
            # pl.figure(1)
            # pl.plot(datanorm_m,'b')
            # pl.plot(self.PCALIB[molecule]['PCA']['norm_interp'][mask,0],'r')
            # pl.plot(pc2,'g')
            # # pl.plot(pc2_inv,'y')
            # #
            # # pl.figure(2)
            # # pl.hist(sqrt((datanorm_m-pc2)**2)/len(datanorm[mask]),100)
            # # # # pl.scatter(self.PCALIB[molecule]['PCA']['norm_interp'][mask,1],(sqrt((datanorm[mask]-self.PCALIB[molecule]['PCA']['norm_interp'][mask,1]))**2))
            # pl.show()


    def rank_molecules(self):

        molkeys = self.PCALIB.keys()
        distance = []
        for molecule in molkeys:
            distance.append(self.PCALIB[molecule]['preselect']['euclid_dist'])

        idx = np.argsort(distance)

        diff = 0
        diffidx = 0
        sortdist = np.asarray(distance)[idx]
        for i in range(len(sortdist)-1):
            if (sortdist[i+1]-sortdist[i]) > diff:
                diff = (sortdist[i+1]-sortdist[i])
                diffidx = i

        self.mol_rank = np.asarray(molkeys)[idx]
        self.mol_dist = sortdist
        if diffidx < 2:
            self.mol_idx = 3
        else:
            self.mol_idx  = diffidx

        # print ''
        # print distance
        print np.asarray(distance)[idx]
        # print np.asarray(molkeys)
        # print ''
        print np.asarray(molkeys)[idx]
        # print diffidx, diff
        #
        # pl.figure(3)
        # pl.plot(np.asarray(distance)[idx])
        # show()


    def calc_astroparams(self):
    #calculating planetary parameters from stellar and orbital parameters

        #calculating mean planetary surface temperature
        self.Tplanet = self.params.star_temp * sqrt(self.params.star_radius / (
            2. * self.params.planet_sma)) * (1 - self.params.planet_albedo) ** (1. / 4.)
        # self.Tplanet = 1400


    def update_params(self):
    #updates the parameter object with perselector derived values and returns
    #updated copy to main code

        newparams = deepcopy(self.params)

        #setting useATM_file to False
        newparams.in_use_ATMfile = False

        #setting planetary temperature
        newparams.planet_temp = self.Tplanet

        #setting number of gases/molecules
        newparams.tp_num_gas = int(self.mol_idx+1)

        #setting molecules list
        newparams.planet_molec = self.mol_rank[self.mol_idx]

        #setting new abs_files path
        newparams.in_abs_path = self.params.pre_cross_path

        #determining correct abs files to be read in
        #reading available cross section lists in PATH
        globlist = glob.glob(self.params.pre_cross_path+'*.abs')

        #determining the right abs file for correct Tplanet
        #this needs to be changed when we want several temperatures
        absfilelist = []
        for molecule in (self.mol_rank[:self.mol_idx+1]):
            temp = self.PCALIB[molecule]['temps']
            next_temp = find_nearest(temp,self.Tplanet)[0]

            for FILE in globlist:
                fname = string.rsplit(FILE,'/',1)[1] #splitting the name
                splitname = string.split(fname,'_',3)

                if splitname[0] == molecule:
                    if float(splitname[2][:-1]) == next_temp:
                        absfilelist.append(fname)

        #setting new list of abs files to be read
        newparams.in_abs_files = absfilelist

        return newparams
