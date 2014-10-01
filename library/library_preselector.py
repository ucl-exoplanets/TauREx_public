import numpy as np
import pylab as pl
import sklearn.decomposition as sk
import glob,string,gzip,os
import cPickle as pickle

from classes.transmission import *
from classes.emission import *
from classes.profile import *
from classes.data import *



def convert2microns(PATH, upcut=25):
#Function converting ExoMol cross section files in dir:PATH from wavenumbers to microns and sorting
#with ascending wavelength
#output: .abs files

    FILES = glob.glob(PATH)

    for f in FILES:
        if os.path.isfile(f[:-3]+'abs'):
            pass
        else:
            print 'converting: ', f
            tmp = np.loadtxt(f,dtype=np.float32)[1:,:]
            tmp[:,0] = 10000.0/tmp[:,0]
            idx = np.argsort(tmp[:,0],axis=-1)
            tmp2 = tmp[idx,:][np.where(tmp[idx,0] < upcut)]
            tmp2 = tmp2.astype(np.float32,copy=False)
            np.savetxt(f[:-3]+'abs',tmp2,fmt="%.6e,%.8e")


def generate_spectra_lib(PARAMS,PATH,OUTPATH,MIXING=[1e-6,1e-5,1e-4,1e-3,1e-2]):
    #Generates transmission spectra from cross section files for given mixing ratios
    #output: 2 column ascii files (.spec)
    
    #cleaning speclib folder before generating spectra
    oldfiles = glob.glob(OUTPATH+'*.spec')
    for i in oldfiles:
        os.remove(i)

    #initiating objects needed
    dataob_pca = data(PARAMS)
    dataob_pca.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)

    profileob_pca = profile(PARAMS, dataob_pca)
    transob_pca = transmission(PARAMS, dataob_pca)

    #reading available cross section lists in PATH
    globlist = glob.glob(PATH+'*.abs')

    if not os.path.isdir(OUTPATH):
        os.mkdir(OUTPATH)
    # else:


    for FILE in globlist:
        print FILE
        fname = string.rsplit(FILE,'/',1)[1] #splitting the name
        temp  = float(string.rsplit(fname,'_',2)[1][:-1]) #getting temperature from file name

        # print fname
        if PARAMS.pre_restrict_temp: #imposing temperature range restrictions
            if temp < float(PARAMS.pre_temp_range[0]) or temp > float(PARAMS.pre_temp_range[1]): 
                run = False
            else: 
                run = True
                
        if run:
            for mix in MIXING:
                X_in   = zeros((1,int(profileob_pca.nlayers))) #setting up mixing ratio array
                X_in  += mix #setting mixing ratio
                rho_in = profileob_pca.get_rho(T=temp) #calculating T-P profile
            
                if PARAMS.in_include_cia or PARAMS.in_include_cld:
                    dataob_pca.set_ABSfile(path=PATH,filelist=[fname],interpolate=True) #reading in cross section file
                else:
                    dataob_pca.set_ABSfile(path=PATH,filelist=[fname],interpolate=False) #reading in cross section file
                        # dataob_pca.specgrid = dataob_pca.wavegrid
                        # dataob_pca.nspecgrid = dataob_pca.nwave
                transob_pca.reset(dataob_pca) #resets transob to reflect changes in dataob
    #                     pl.figure(21)
    #                     pl.plot(dataob_pca.specgrid,dataob_pca.sigma_array[0])
    #                     pl.show()
    #         
                        #manually setting mixing ratio and T-P profile
                MODEL = transob_pca.cpath_integral(rho=rho_in,X=X_in) #computing transmission
#                 pl.figure(20)
#                 pl.plot(MODEL)
#                 pl.show()
    #               exit()
                np.savetxt(OUTPATH+fname[:-4]+'_'+str(mix)+'d.spec',np.column_stack((dataob_pca.specgrid,MODEL)))


def find_nearest(arr, value):
    # find nearest value in array
    arr = np.array(arr)
    idx = (abs(arr-value)).argmin()
    return [arr[idx], idx]


def PCA(DATA0,PCnum):
#small wrapper doing PCA using sk module
    DATA = np.transpose(DATA0)
    # DATA = DATA0

    u,S,PC = np.linalg.svd(DATA)
    feature = np.transpose(PC)
    print 'no. of PCs: ', len(PC[:,0])
    feature = feature[:PCnum,:]
    return feature, PC



def generate_PCA_library(PARAMS,PATH,OUTPATH=False,comp_num=2):
        
    FILES = glob.glob(PATH)

    filename = FILES[0].rpartition('/')[-1]
    filenamesplit = filename.split('_')
    molname = filenamesplit[0]

    molnamelist = []
    molnamelist.append(molname)
    molnumlist = []
    mollenlist = []
    for f in FILES:
        if (f.rpartition('/')[-1]).split('_')[0] != molnamelist[-1]:
            molnamelist.append((f.rpartition('/')[-1]).split('_')[0])


    for i in range(len(molnamelist)):
        count = 0
        lencount = 0
        idx = 0
        for f in FILES:
            if (f.rpartition('/')[-1]).split('_')[0] == molnamelist[i]:
                count += 1
                if idx == 0:
                    lencount = np.shape(np.loadtxt(f))[0]
                    idx = 1

        mollenlist.append(lencount)
        molnumlist.append(count)

    OUTdic = {}
    for i in range(len(molnamelist)):
        DATA = np.zeros((mollenlist[i],molnumlist[i]))

        templist = []
        j=0
        for f in FILES:
            if (f.rpartition('/')[-1]).split('_')[0] == molnamelist[i]:
                tmp = np.loadtxt(f)
                temp  = float(string.rsplit(f,'_',3)[1][:-1]) #getting temperature from file name
                
                if PARAMS.pre_restrict_temp: #imposing temperature range restrictions
                    if temp > float(PARAMS.pre_temp_range[0]) and temp < float(PARAMS.pre_temp_range[1]): 
                        print f
                        DATA[:,j] = tmp[:,1]
                        if temp not in templist:
                            templist.append(temp)
                else:
                    DATA[:,j] = tmp[:,1]
                    if temp not in templist:
                        templist.append(temp)

                if j == 0:
                    wavegrid = tmp[:,0]
                j += 1

        pca = sk.RandomizedPCA(n_components=comp_num,whiten=False)
        
        #checking and replacing negative numbers 
        DATA[DATA < 0] = 0
        
#         pl.figure(1)
#         pl.plot(DATA)
#         pl.show()
        
        pca.fit(np.transpose(DATA))

        meanspec = np.mean(DATA,axis=1)

        normPCA = np.zeros((mollenlist[i],comp_num))
        for jj in range(comp_num):
            normPCA[:,jj] = (pca.components_[jj]-np.min(pca.components_[jj])) /np.max((pca.components_[jj]-np.min(pca.components_[jj])))
        
#         print pca.explained_variance_ratio_ 
#         pl.figure(1)
#         pl.plot(pca.components_[0],'b')
#         pl.figure(2)
#         pl.plot(pca.components_[1],'r')
#         pl.show()
        
        
#         pl.figure(2)
#         pl.plot(normPCA[:,0],'b')
#         pl.plot(normPCA[:,1],'r')
        

        #add to OUTdic
        OUTdic[molnamelist[i]] = {}
        OUTdic[molnamelist[i]]['length']   = int(mollenlist[i])
        OUTdic[molnamelist[i]]['numbers']  = int(molnumlist[i])
        OUTdic[molnamelist[i]]['wavegrid'] = wavegrid
        OUTdic[molnamelist[i]]['temps']    = sorted(templist)
        OUTdic[molnamelist[i]]['data']     = DATA
        OUTdic[molnamelist[i]]['meanspec'] = meanspec
        OUTdic[molnamelist[i]]['PCA']      = {}
        OUTdic[molnamelist[i]]['PCA']['full'] = pca.components_
        OUTdic[molnamelist[i]]['PCA']['norm'] = normPCA
        # OUTdic[molnamelist[i]]['filename'] =
        # OUTdic[molnamelist[i]]['path']     = PATH



    if OUTPATH != False:
        with gzip.GzipFile(OUTPATH+'spec_pcalib.pkl.zip','wb') as outhandle:
            pickle.dump(OUTdic,outhandle)

    return OUTdic