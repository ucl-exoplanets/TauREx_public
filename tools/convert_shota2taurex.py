#! /usr/bin/python 

##############
# Small script that converts the Input files by Shota Notsu to generic TauREx format
##############

import numpy as np
import argparse 


class convert_shota(object):
    def __init__(self,options): 
        #initialising object 
        
        self.path = options.inpath
        self.outpath = options.outpath
        
        # np.loadtxt(path)
        
        #Lists of possible molecules and corresponding moleuclar weights. Just add to the lists if molecules are not available. 
        self.molnames = ['C_1D','H','N','O','O_1D','O_1S', 'CO', 'H2','HO','N2', 'NO', 'O2', 'O2_D', 'O3', 'CH4', 
                    'CO2', 'H2O', 'HO2', 'N2O', 'NO2', 'H2O2', 'HNO3', 'CH2O2', 'HCOOH', 'CH3ONO', 'e-', 
                    'H+', 'O+', 'NO+','O2+', 'C','HN','CNC','H2N','H3N','C+','C-','N+','O-','CO+','HO+','N2+','CHO+',
                    'CH3','CHO','HCN','HNO','NO3','C2H2','C2H6','CH2O','HNO2','N2O3','CH3O2','CH3OH','CH4O2','H3O+','HE','NH3']
        self.molweights = [14,1,14,16,18,18,28,2,17,28,30,32,34,48,16,44,18,33,44,46,34,63,46,46,61,0,1,16,30,32,12,15,38,16,17,12,12,14,16,28,17,28,29,
                      15,29,27,31,62,26,28,30,48,76,47,32,52,19,4,17]
    
    
        #loading data
        self.load_data()
        #assigning molecules 
        self.assign_molecules()
        #saving data 
        self.save_data()
        
    
    
    def load_data(self):
        #reading the molecule list in input file
        with open(self.path, 'r') as f:
                        for i in range(4):
                            f.readline()
                        self.molecules_raw = [x.upper() for x in f.readline().split()]
        
        #loading rest of data 
        self.data = np.loadtxt(self.path,skiprows=7)[::-1]

    def assign_molecules(self):
        #assigning molecular weights to molecules
        self.mollist = []
        self.weightlist = []
        for molecule in self.molecules_raw:
            mol = molecule.split('_')[0]
            self.mollist.append(mol)
            if mol in self.molnames:
                self.weightlist.append(self.molweights[self.molnames.index(mol)])
            else: 
                print('NOT FOUND: ',mol)
        
    
        print('Molecules assigned: ',self.mollist)
        print('Molecular weights: ',self.weightlist)
    
    def save_data(self):
        with open(self.outpath,'wb') as outfile:
            for mol in self.mollist:
                outfile.write(mol+' ')
            outfile.write('\n')
            for weight in self.weightlist:
                outfile.write(np.str(weight)+' ')
            outfile.write('\n')
            np.savetxt(outfile,self.data)
    

if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--infile',
                      dest='inpath',
                      default='',
                       help='Input file'
                      )
    parser.add_argument('--outfile',
                        dest='outpath',
                        default='abundance.dat',
                        help='output filename')

    
    #parsing command line options
    options = parser.parse_args()
    
    #loading main object 
    convert_shota(options)
    

