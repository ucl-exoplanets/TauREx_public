#! /usr/bin/python

################################################
#convertwavenumbers.py
#
# Small script converting 2 column ExoMol cross-sections
# from wavenumbers to microns
#
# Input: - dir path of cross section folder
#
#
# Output: - .abs files (same than input but with microns in first column)
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Feb 2014
#
################################################

import numpy as np
import glob
import os


CWD = os.getcwd()
PATH = CWD+'/crosssections/*'

#cutoffs in microns otherwise gets unnecessarily big
lowcut = 0;
upcut = 25;


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


