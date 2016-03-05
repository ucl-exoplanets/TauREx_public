#! /usr/bin/python

#small script running analayse_from_traces on many folders (each folder is taurex output, i.e. stage_0/1)
#if more than one folder is found, the script will divide the tasks onto multiple cores

import os, glob, optparse
import numpy as np
import multiprocessing as mp

parser = optparse.OptionParser()
parser.add_option('-d', '--dir',
                  dest="path",
                  default="Output",
                  )
parser.add_option('-n', '--nprocess',
                  dest="max_nprocesses",
                  default=1,
                  )
parser.add_option('-p', '--parfile',
                  dest='parfile',
                  default=None)


options, remainder = parser.parse_args()
# options.path = '/Volumes/DATA_CALIMERO/ingo/taurex2_paper/revised'



folder_list = glob.glob(options.path+'/*')
N_folder = len(folder_list)

N_cpu = mp.cpu_count()
print 'No. of Folders to process: ', N_folder

if N_folder >= N_cpu:
    N_processes = N_cpu
else:
    N_processes = N_folder

if int(N_processes) > int(options.max_nprocesses):
    N_processes = options.max_nprocesses

print 'No. of CPUs available: ', N_cpu
print 'No. of CPUs used: ',N_processes

print 'Starting plotting... '

# print N_processes
print folder_list

def plot_folder(folder):
    if os.path.isdir(folder):
        print 'Plotting: ', folder
        if options.parfile is None:
            parfile = glob.glob(folder+'/*.par')[0]
        else:
            parfile = options.parfile
        os.system('python analyse_solutions_from_traces.py -p '+parfile+' -d '+folder)
        os.system('python taurex_mlehess.py -p '+parfile+' -d '+folder)

pool = mp.Pool(processes=int(N_processes))           #setting number of cores on which to run
pool_result = pool.map(plot_folder,folder_list) #runnning the stuff
pool.close()                                    #closing the pool 
# pool.join()                                     #closing the pool.



##########old loop
# for folder in folder_list:
#     if os.path.isdir(folder):
#         print 'Plotting: ', folder
#         parfile = glob.glob(folder+'/stage_0/*.par')[0]
#         os.system('python analyse_solutions_from_traces.py -p '+parfile+' -d '+folder)
# #     exit()