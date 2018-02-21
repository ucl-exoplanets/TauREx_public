#! /usr/bin/python 
'''
Script converting Venot abundance and pressure files to Tau-REx input files. 
Tau-REx input for abundances: 

Molecule names
Pressure (Pa), Abundances... 

Tau-REx input for TP-profiles:

Pressure (Pa), Temperature 
'''

import numpy as np 
import pylab as pl
import argparse
import os


class convert_venot(object):
    def __init__(self,options): 
        
        self.options = options
        self.type = options.type 
        self.filename = self.options.input
        if self.type == 'abundance':
            self.load_abundance()
            self.convert_abundance()
        elif self.type == 'temperature':
            self.load_tp()
            self.convert_tp()
        
        self.save_file()
        
        
    def load_abundance(self):
        '''
        function loading venot abundance file 
        '''
        with open(self.filename, 'r') as f:
            first_line  = f.readline()
            second_line = f.readline()
            self.molecules = ' '.join(first_line.split()).split(' ')
            self.mmw       = ' '.join(second_line.split()).split(' ')

        # read mixing ratio profiles (from third row)
        self.mixratio_profiles_file = np.loadtxt(self.filename, skiprows=2)
     
        print self.molecules
        print self.mmw
        exit()
        
    def load_tp(self):
        '''
        function loading venot temperature-pressure file
        '''
    def convert_abundance(self):
        '''
        function converts venot file to taurex input file. 
        Venot: 
        Molecules
        Molecular weights
        1st column, Pressure (Pa), Abundances 
        
        Taurex:
        Molecules
        Pressure (Pa), Abundances
        '''
        
        ble = 1
        
    def convert_tp(self):
        '''
        function converts venot TP profile file to taurex input file. 
        Venot: 
        1st colum, Pressure (Pa), Temperature (K), 4th column, 5th column 
        
        Taurex:
        Pressure (Pa), Temperature (K)
        '''
        
        ble =1 
        
    def save_file(self):
        '''
        function saving file given type 
        '''

        
if __name__ == '__main__':

    #loading parameter file parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',
                      dest='input',
                      default='',
                       help='Input file'
                      )
    parser.add_argument('--type',
                        dest='type',
                        default='abundance',
                        help='Type of file: abundance / temperature')
    parser.add_argument('--outfile',
                        dest='outfile',
                        default='abundance.dat',
                        help='output filename')
    parser.add_argument('--outdir',
                        dest='outdir',
                        default='',
                        help='output directory')
    
    #parsing command line options
    options = parser.parse_args()
    
    
    #initialising binning object 
    convert_ob = convert_venot(options)
    
    
    #saving converted file
#     print 'saving file to: {}'.format(os.path.join(options.out_dir,options.out_file))
#     convert_ob.save_file(options.out_dir, options.out_file)
        
 