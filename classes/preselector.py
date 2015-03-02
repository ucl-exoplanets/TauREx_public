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
from base import base
import numpy as np
import pylab as pl
import gzip,os, copy
import cPickle as pickle
# from copy import deepcopy
import library_preselector as lib_pre
# from library_preselector import *
import scipy.stats.stats as st
import logging

class preselector(base):

    def __init__(self, model_object, data=None, params=None):

        logging.info('Initialising Marple object')


        if params:
            self.params = params
        else:
            self.params = model_object.params

        if data:
            self.data = data
        else:
            self.data = model_object.data

        self.model_object = model_object

        #setting emission/transmission model
        self.set_model(model_object)

        #reading in spectrum data to be fitted
        self.spectrum = self.data.spectrum
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
        self.generate_mask(thres=self.params.pre_mask_thres)
        #computing euclidian distance and correlation coefficients between library and data
        self.correlate()
        #rank molecules according to their lowest euclidian distance to data
        #also estimate number of likely molecules in data
        self.rank_molecules()
        #selecting molecules 
        self.select_molecules()
        #calculating some planetary parameters
        self.calc_astroparams()
        
        #small hack for paper excluding water
#         print 'moleselected ',self.molselected
#         print 'numgas ', self.numgas
#         
#         self.molselected= np.delete(self.molselected, 2)
#         self.numgas -= 1
#         print 'newmolselected ', self.molselected
#         print 'newnumgas ', self.numgas
        
#         exit()

        # print self.mol_rank
        # print self.mol_idx
        
        #saving new PCA library
#         with gzip.GzipFile(self.params.pre_pca_path+'spec_pcalib2.pkl.zip','wb') as outhandle:
#             pickle.dump(self.PCALIB,outhandle)



    def run_preprocess(self,generateSpectra=None,generatePCA=None):
        #function running pre-processing steps
        #convertLinelist will convert ExoMol cross sections from wavenumbers to microns
        #generateSpectra will generate library of transmission spectra from linelists
        #generatePCA will run PCA on spectral library and compile final dictionary

        #getting booleans from parameter file if not set manually
        if generateSpectra is None:
            generateSpectra = self.params.pre_gen_speclib
        if generatePCA is None:
            generatePCA = self.params.pre_gen_pca

        #doing the pre-processing
#             print '1 done'
        if generateSpectra:
            lib_pre.generate_spectra_lib(self.params,self.params.in_abs_path,self.params.pre_speclib_path,
                                 MODEL=self.model_object,MIXING=self.params.pre_mixing_ratios)
            self.params.console.setLevel(20)
#             print '2 done'
        if generatePCA:
            self.params.console.setLevel(30)
            lib_pre.generate_PCA_library(self.params,self.params.pre_speclib_path+'*',self.params.pre_pca_path)
            self.params.console.setLevel(20)
#             print '3 done'


    def load_library(self,PATH=None):
        #function loading spec_pcalib.pkl.zip into memory
        if PATH is None:
            PATH = self.params.pre_pca_path

        try:
            with gzip.open(PATH+'spec_pcalib.pkl.zip',mode='rb') as filehandle:
                self.PCALIB = pickle.load(filehandle)
        except IOError:
            logging.warning('WARNING: cannot find spec_pcalib.pkl.zip in: ',PATH)
            redo = raw_input('Try to generate library from scratch? (y/n) [n]: ')
            if redo == 'N' or redo == 'n' or redo == 'no' or redo == '':
                exit()
            elif redo == 'Y' or redo == 'y' or redo == 'yes':
                logging.info('generating library...')
                self.run_preprocess(convertLinelist=True,generateSpectra=True,generatePCA=True)
                logging.info('loading library...')
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
                
                #applying a second normalisation to be within 0-1 of self.wavegrid range
                interpc[:,i] = (interpc[:,i] - np.min(interpc[:,i])) / np.max((interpc[:,i] - np.min(interpc[:,i])))
                
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
                self.PCALIB[molecule]['PCA']['interp_mask'] = pc1_interp > self.params.pre_mask_thres
            except(KeyError):
                pass


            # print mask
#             print molecule
#             mask = pc1 > thres
#             pc2 = self.PCALIB[molecule]['PCA']['norm'][:,1]
#             xaxis = np.zeros((len(pc1),1))
#             for j in range(len(pc1)):
#                 xaxis[j] = j
#             # # # print pc1
#             # # # exit()
#             pl.figure(1)
#             plot(self.PCALIB[molecule]['wavegrid'],self.PCALIB[molecule]['PCA']['norm'])
#             pl.plot(xaxis,pc1,'b')
#             pl.plot(xaxis,pc2,'g')
#             pl.plot(xaxis[mask],pc1[mask],'r')
#             # #
#             pl.figure(2)


#             pl.figure(1)
#             self.PCALIB[molecule]['PCA']['norm_interp'][:,0]
#             pl.figure(2)
#             pl.plot(self.PCALIB[molecule]['PCA']['interp_mask'],'r')
# 
#             pl.show()
            # # exit()

    def correlate(self):
        #applies pre-determined masks and correlates the data to the second
        #principal component of the spectral library

        #normalise data
        data = self.spectrum[:,1]
        datanorm = data - np.min(data)
        datanorm /= np.max(datanorm)
        errnorm = self.spectrum[:,2] * (datanorm / (data-np.min(data)))
#         errnorm = errnorm[:-1]


        molkeys = self.PCALIB.keys()
        for molecule in molkeys:
            #adding layer to PCALIB
            self.PCALIB[molecule]['preselect'] = {}

            #loading mask
            mask = self.PCALIB[molecule]['PCA']['interp_mask']

            datanorm_m = datanorm[mask]-np.mean(datanorm[mask])

            pc1 = self.PCALIB[molecule]['PCA']['norm_interp'][mask,0] - np.mean(self.PCALIB[molecule]['PCA']['norm_interp'][mask,0])
            pc2 = self.PCALIB[molecule]['PCA']['norm_interp'][mask,1] - np.mean(self.PCALIB[molecule]['PCA']['norm_interp'][mask,1])

            pc1_inv = (-1.0*(pc1-np.mean(pc1)))#+np.mean(pc1)
            pc2_inv = (-1.0*(pc2-np.mean(pc2)))#+np.mean(pc2)


            # calculating the euclidian distance between PC1/PC2 and normalised data 
            try:
                dist1 = np.abs((datanorm_m-pc1))
                eucdist_pc1 = np.sum(dist1) /  float(len(datanorm)) #len(datanorm[mask])* (float(len(datanorm))/float(len(datanorm[mask])))
            except ZeroDivisionError:
                eucdist_pc1 = 1000
            try:
                dist2 = np.abs((datanorm_m-pc2))
                eucdist_pc2 = np.sum(dist2)/  float(len(datanorm)) # len(datanorm[mask])* (float(len(datanorm))/float(len(datanorm[mask])))
            except ZeroDivisionError:
                eucdist_pc2 = 1000
                
            try:
                dist3 = np.abs((datanorm_m-pc1_inv))
                eucdist_pc1_inv = np.sum(dist3)/  float(len(datanorm)) #len(datanorm[mask])* (float(len(datanorm))/float(len(datanorm[mask])))
            except ZeroDivisionError:
                eucdist_pc1_inv = 1000
            try:
                dist4 = np.abs((datanorm_m-pc2_inv))
                eucdist_pc2_inv = np.sum(dist4)/  float(len(datanorm)) #len(datanorm[mask])* (float(len(datanorm))/float(len(datanorm[mask])))
            except ZeroDivisionError:
                eucdist_pc2_inv = 1000
                
            #checking for broken numbers
            if np.isnan(eucdist_pc1):
                eucdist_pc1 = 1000
            if np.isnan(eucdist_pc1_inv):
                eucdist_pc1_inv = 1000
            if np.isnan(eucdist_pc2):
                eucdist_pc2 = 1000
            if np.isnan(eucdist_pc2_inv):
                eucdist_pc2_inv = 1000

#             eucdist = [eucdist_pc1,eucdist_pc2] #include both principal components for correlation analysis
            eucdist = [eucdist_pc1]  #only include the first PC
#             eucdist_inv = [eucdist_pc1_inv, eucdist_pc2_inv] #include both PCs for correlation, flipped upside down
            eucdist_inv = [eucdist_pc1_inv] #only include the first PC
            
#             logging.info('Molecule: %s, eucdist: %d, inv_eucdist: %d' % (molecule,eucdist[0],eucdist_inv[0]))
#             print 'molecule ', molecule
#             print 'eucdist ',eucdist
#             print 'eucdist_inv ', eucdist_inv
            
            #calculating pearson correlation coefficients between PCs and data 
            corrcoeff_pc1 = st.pearsonr(datanorm_m,pc1)
            corrcoeff_pc2 = st.pearsonr(datanorm_m,pc2)
            


            self.PCALIB[molecule]['preselect']['pearson'] = corrcoeff_pc1
            self.PCALIB[molecule]['preselect']['euclid_dist'] = min(eucdist)
            
            #@todo the inverse correlation is now disabled. doesnt make too much sense for transmission but may need implementation for emission.
#             if corrcoeff_pc1[0] <0.0:
# #             if min(eucdist_inv) < min(eucdist):
#                 self.PCALIB[molecule]['preselect']['euclid_dist'] = min(eucdist_inv)
# 
#             else:
#                 self.PCALIB[molecule]['preselect']['euclid_dist'] = min(eucdist)
                


            logging.info('Molecule: {:10s}, corrcoeff: {:4g}, eucdist: {:4g}, inv_eucdist: {:4g}'.format(molecule,corrcoeff_pc1[0], np.min(eucdist),np.min(eucdist_inv)))
            
#             print molecule,': ',corrcoeff_pc1, '... ',eucdist,'... ',eucdist_inv
            xnums = np.arange(len(datanorm_m))
            xnums_pc1 = np.arange(len(pc1))
  
#             pl.figure(1)
#             pl.plot(xnums,datanorm_m,'b')
# #             pl.errorbar(xnums,datanorm_m,errnorm[:-1])
#             pl.plot(xnums_pc1,pc1,'r')
#             pl.plot(xnums_pc1,pc2,'g')
            
#             
#             pl.figure(2)
#             pl.plot(self.PCALIB[molecule]['wavegrid'],self.PCALIB[molecule]['PCA']['norm'][:,0],'r')
#             pl.plot(self.PCALIB[molecule]['wavegrid'],self.PCALIB[molecule]['PCA']['norm'][:,1],'g')
#             pl.xlim([self.wavegrid[0], self.wavegrid[-1]])
# # #             # pl.plot(pc2_inv,'y')
# # #             #
# # #             pl.figure(2)
# # #             pl.hist(sqrt((datanorm_m-pc2)**2)/len(datanorm[mask]),100)
# #             # # # pl.scatter(self.PCALIB[molecule]['PCA']['norm_interp'][mask,1],(sqrt((datanorm[mask]-self.PCALIB[molecule]['PCA']['norm_interp'][mask,1]))**2))
#             pl.show()


    def rank_molecules(self):

        molkeys = self.PCALIB.keys()
        distance = []
        for molecule in molkeys:
            print molecule, self.PCALIB[molecule]['preselect']['euclid_dist']
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
        print self.mol_rank
        self.mol_dist = sortdist
        if diffidx < 2:
            self.mol_idx = 4
        else:
            self.mol_idx  = diffidx
        self.mol_idx = diffidx
        
#         self.mol_idx = 5

        # print ''
        # print distance
#         print np.asarray(distance)[idx]
        # print np.asarray(molkeys)
        # print ''
#         print np.asarray(molkeys)[idx]
        # print diffidx, diff
        #
#         pl.figure(3)
#         pl.plot(np.asarray(distance)[idx])
#         pl.show()


    def calc_astroparams(self):
    #calculating planetary parameters from stellar and orbital parameters

        #calculating mean planetary surface temperature
        self.Tplanet = self.params.star_temp * np.sqrt(self.params.star_radius / (
            2. * self.params.planet_sma)) * (1 - self.params.planet_albedo) ** (1. / 4.)
        # self.Tplanet = 1400


    def select_molecules(self):
    #select molecules given previously determined ranking and inclusions from parameter file
        
        self.molselected = self.mol_rank[:self.mol_idx]
        tmpmol = self.molselected
#         self.numgas      = int(self.mol_idx+1)

        #if set in parameters, molecules are added to selection if not already present
        if self.params.pre_mol_force_bool:
            if np.ndim(self.params.pre_mol_force) == 0:
                self.params.pre_mol_force = np.array([self.params.pre_mol_force])

            for mol in self.params.pre_mol_force:
                if mol.strip() not in self.molselected:
                    tmpmol = np.insert(tmpmol,0,mol.strip())
#                     self.numgas += 1  
                 
            self.molselected = tmpmol
            
        self.numgas = int(len(self.molselected))


            
        

    def update_params(self):
    #updates the parameter object with perselector derived values and returns
    #updated copy to main code

        # @todo looks like only important parameter that changes is params.planet_molec

#         logging.info('Update parameters object: %s' % self.molselected)

        
#         newparams = copy.deppcopy(self.params) #@todo due to the logging in the parameter class we cannot copy it anymore. 
        newparams = self.params #@todo setting parameters directly in original instance. may not matter or be a very bad idea. Ideas?

        #setting useATM_file to False
        newparams.in_use_ATMfile = False # @todo really needed? - IPW: yes i think so 

        #setting planetary temperature
#         newparams.planet_temp = self.Tplanet

        #setting number of gases/molecules
        newparams.ngas = self.numgas # @todo deprecated!

        #setting molecules list
        newparams.planet_molec = self.molselected

        #setting new abs_files path
        newparams.in_abs_path = self.params.in_abs_path # @todo why? Does it change? -IPW: it used to but not any more so yes... useless

#         #determining correct abs files to be read in
#         #reading available cross section lists in PATH
#         globlist = glob.glob(self.params.in_abs_path+'*.abs')
# 
#         #determining the right abs file for correct Tplanet
#         #this needs to be changed when we want several temperatures
#         absfilelist = []
#         for molecule in (self.molselected):
#             temp = self.PCALIB[molecule]['temps']
#             next_temp = find_nearest(temp,self.Tplanet)[0]
# 
#             for FILE in globlist:
#                 fname = string.rsplit(FILE,'/',1)[1] #splitting the name
#                 splitname = string.split(fname,'_',3)
# 
#                 if splitname[0] == molecule:
#                     if float(splitname[2][:-1]) == next_temp:
#                         absfilelist.append(fname)
# 
#         #setting new list of abs files to be read
#         newparams.in_abs_files = absfilelist

        return newparams
