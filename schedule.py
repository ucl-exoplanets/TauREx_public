#! /usr/bin/python

import numpy as np
import os, glob 


outdir = '/Volumes/DATA/Ingo/taurex/ariel/'

# list = ['r300_1e-5','r300_5e-5','r300_1e-4','r300_5e-4','r30_5e-4','r30_1e-4']
# list = ['r200_1e-5','r200_5e-5','r200_1e-4','r200_5e-4','r100_5e-5','r50_5e-5','r30_5e-5','r300_5e-5','r100_1e-5','r50_1e-5','r30_1e-5','r300_1e-5']

# list = ['r50_snr_7','r50_snr_10','r50_snr_20','r100_snr_7','r100_snr_10','r100_snr_20','r200_snr_7','r200_snr_10','r200_snr_20']
# 
# for name in list:
#     print '------------------------------------------'
#     print name
#     print ''
#     os.system('mpirun -np 6 python exonest.py -p Parfiles/ariel/exonest_hotjupiter_ariel_'+name+'.par')
#     os.system('cp Output/model_fit.pdf chains/')
#     os.system('cp Parfiles/exonest_hotjupiter_'+name+'.par chains/')
#     os.system('cp -r chains/ '+outdir+'chains_ariel_hotjupiter_'+name+'')

list = ['r20_snr_7','r20_snr_10','r50_snr_7','r50_snr_10','r100_snr_7','r100_snr_10']

for name in list:
    print '------------------------------------------'
    print name
    print ''
    os.system('mpirun -np 4 python exonest.py -p Parfiles/ariel/exonest_warmneptune_ariel_'+name+'.par')
#     os.system('cp Output/model_fit.pdf chains/')
    os.system('cp Parfiles/ariel/exonest_warmneptune_'+name+'.par chains/')
    os.system('cp -r chains/ '+outdir+'chains_ariel_warmneptune_'+name+'')