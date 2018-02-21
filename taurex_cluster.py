#! /usr/bin/python -W ignore

###########################################################
# taurex_cluster
#
# submits and runs multiple instances of taurex on cobweb/legion.
# 
# Input: 
#      -d command flag: ascii dictionary name required by cluster class
#      all other inputs defined by dictionary. See create_cluster_dict
#
# Output: 
#      standard Taurex output in variable output folders definable by dictionary
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, May 2015
#
###########################################################


import  sys, os, optparse, time
import numpy as np

sys.path.append('./classes')
# import cluster
from cluster import cluster

#loading parameter file parser
parser = optparse.OptionParser()
parser.add_option('-d', '--dictname',
                  dest="dictname",
                  default="run_taurex.dict",
                  action="store"
                  )
parser.add_option('-i', '--datadir',
                  dest="datadir",
                  default="pwd",
                  action="store"
                  )
                  

options, remainder = parser.parse_args()

if options.datadir is "pwd":
    options.datadir = os.getcwd()


#loading cluster object
c = cluster(options=options)

#reading dictionary
c.read_dict()

#writing run script and submitting to queue
for IDs in c.IDs:
    scriptname = c.generate_script(IDs)
#     os.system('qsub {}'.format(scriptname))  #script submit command for cobweb/legion
#     os.system('rm {}'.format(scriptname))
#     time.sleep(5)

#   the equivalent if to be run on pingu
#     os.system('mpirun -np 12 python taurex.py -p Parfiles/taurex_emission_wasp76.par -c run_taurex.dict -i '+str(IDs))
#     os.system('mkdir -p /Volumes/DATA_PINGU/ingo/cobweb/tests/pingu_wasp76paper_10000/'+str(IDs))
#     os.system('cp -rf /Users/ingo/repos/taurex/Output/* /Volumes/DATA_PINGU/ingo/cobweb/tests/pingu_wasp76paper_10000/'+str(IDs))


