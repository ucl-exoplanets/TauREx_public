#! /usr/bin/python

import numpy as np
import os, glob 


outdir = '/Volumes/DATA/Ingo/taurex/'

# list = ['r300_1e-5','r300_5e-5','r300_1e-4','r300_5e-4','r30_5e-4','r30_1e-4']
list = ['r200_1e-5','r200_5e-5','r200_1e-4','r200_5e-4','r100_5e-5','r50_5e-5','r30_5e-5','r300_5e-5','r100_1e-5','r50_1e-5','r30_1e-5','r300_1e-5']


for name in list:
    print '------------------------------------------'
    print name
    print ''
    os.system('mpirun -np 8 python exonest.py -p Parfiles/exonest_hotjupiter_'+name+'.par')
    os.system('cp Output/model_fit.pdf chains/')
    os.system('cp Parfiles/exonest_hotjupiter_'+name+'.par chains/')
    os.system('cp -r chains/ '+outdir+'chains_'+name+'_new')


