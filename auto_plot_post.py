#! /usr/bin/python

#small script running analayse_from_traces on many folders (each folder is taurex output, i.e. stage_0/1)

import os, glob
import numpy as np


path = '/Volumes/DATA_CALIMERO/ingo/taurex2_paper/revised'
parfile = 'default.par'

list = glob.glob(path+'/*')

for folder in list:
    if os.path.isdir(folder):
        print 'Plotting: ', folder
        parfile = glob.glob(folder+'/stage_0/*.par')[0]
        os.system('python analyse_solutions_from_traces.py -p '+parfile+' -d '+folder)
#     exit()
