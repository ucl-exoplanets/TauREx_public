#! /usr/bin/python

#small script runnign the plot_chains script on many folders. Used for taurex paper 1 example

import os, glob
import numpy as np


path = '/Users/ingowaldmann/UCLlocal/taurex_paper1/'

list = glob.glob(path+'chains_r*')


for folder in list:
    print 'Plotting: ', folder
#     os.system('../plot_chains.py -n -c False -i '+folder+'/ -o '+folder+'/')
    os.system('./plot_spectrum_from_chains.py -n -p Parfiles/exonest_hotjupiter.par -i '+folder+'/ -o '+folder+'/')
#     exit()