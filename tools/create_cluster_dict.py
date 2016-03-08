#! /usr/bin/python -W ignore

###########################################################
# create_cluster_dict
#
# This is a makeshift script as I'm not sure how to generalise this. 
# It basically generates a dictionary for multiple runs on cobweb/legion. 
# The dictionary is read in by the cluster class, which will generate 
# submit scripts and modify on runtime params object in taurex accordingly.
#
# The exact run setup is likely to be custom built each time. The dictionary has the form:
#
#    dict[ID_number][parameter1][value]    #for parameter object parameters to be modified on runtime
#    dict[ID_number][parameter2][value]
#                                          #general run parameters
#    dict[ID_number][GENERAL][NODES]       #int number of nodes (1-4 cobweb)
#                            [CPUS]        #int number of cores per node (1-24 cobweb) 
#                            [WALLTIME]    #hh:mm:ss
#                            [MEMORY]      #in GB
#                            [PARFILE]     #(e.g) Parfiles/taurex_emission.par 
#                            [OUTPUT_DIR]  #(e.g.) /share/data/taurex/ID_number
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, May 2015
#
###########################################################

import sys,os
import numpy as np
import glob
import string

sys.path.append('./classes')
from cluster import cluster




#parameter values for individual runs 
# alpha_l = [0.0,0.25,0.5,0.75,1.0]  #emission alpha parameter gridded from 0 - 1 in RNUM steps
# alpha_h = [1.0,0.75,0.5,0.25,0.0]

datalist = glob.glob('Input/observations/riemann_hd209/*')

# print datalist

with open('datalist.txt','wb') as ofile:
    for name in datalist:
        ofile.write(string.split(name,'/')[-1]+'\n')

#defining number of iteration per parameter in this case 
RNUM = len(datalist)

#defining general run parameters
GENERAL = {}
GENERAL['NODES']      = 1
GENERAL['CPUS']       = 24
GENERAL['WALLTIME']   = '10:00:00'
GENERAL['MEMORY']     = 100
GENERAL['PARFILE']    = 'Parfiles/riemann_largespec.par'
GENERAL['CLUSTER']    = 'cobweb'
GENERAL['USERNAME']   = 'ucapipw'
GENERAL['DISKMEM']    = 10

if GENERAL['CLUSTER'] is 'legion':
    #legion directories (all must be absolute paths)
    #legion does not require a scratch path, handled by $TMPDIR variable
    GENERAL['DATA_DIR']   = '/home/ucapipw/Scratch/taurex'
    GENERAL['OUTPUT_DIR'] = '/home/ucapipw/Scratch/riemann_output'

elif GENERAL['CLUSTER'] is 'cobweb':
    #cobweb directories (all must be absolute paths)
    GENERAL['DATA_DIR']   = '/share/data/ingo/repos/taurex_cluster'
    GENERAL['SCRATCH_DIR']= '/scratch/ingo/run'
    GENERAL['OUTPUT_DIR'] = '/share/data/ingo/taurex/riemann_hd209'
    


#initialising dictionary
DICT = {}

#filling dict with parameter values per run
ID = 0
for i in range(RNUM): #this may be changed with a more informative ID... e.g. date
    DICT[ID] = {}
    DICT[ID]['GENERAL'] = GENERAL.copy()        #different run setups per run can be provided here
    DICT[ID]['GENERAL']['OUTPUT_DIR'] = GENERAL['OUTPUT_DIR']
    DICT[ID]['in_spectrum_file'] = datalist[i]  #setting parameter for run
#     DICT[ID]['out_path'] = GENERAL['OUTPUT_DIR']+'/{}'.format(ID) #must be provided for legion but not for cobweb
    ID += 1

# ID = 0
# for i in range(RNUM):
#     DICT[ID] = {}
#     DICT[ID]['GENERAL'] = GENERAL.copy()        #different run setups per run can be provided here
#     DICT[ID]['GENERAL']['OUTPUT_DIR'] = GENERAL['OUTPUT_DIR']+'/'+str(ID)
#     DICT[ID]['fit_hybrid_alpha_h'] = alpha_h[i]  #setting parameter for run
#     ID += 1


#writing dictionary to ascii. Use function provided by cluster for compatibility reasons. 
cluster().write_dict(DICT,DICTNAME='run_taurex.dict')

    