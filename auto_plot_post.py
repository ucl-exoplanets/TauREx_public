#! /usr/bin/python

#small script running analayse_from_traces on many folders (each folder is taurex output, i.e. stage_0/1)
#if more than one folder is found, the script will divide the tasks onto multiple cores

import os, glob, optparse, string, functools,sys,logging
import numpy as np
import multiprocessing as mp

import subprocess

    

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
parser.add_option('-c', '--cluster',
                  dest="cluster_dictionary",
                  default="None",
                  type = "string",
                  action="store"         
                  )
# parser.add_option('-i', '--cluster_index',
#                   dest="cluster_procid",
#                   default="None",
#                   type = "string",
#                   action="store"         
#                   )
options, remainder = parser.parse_args()
# options.path = '/Volumes/DATA_CALIMERO/ingo/taurex2_paper/revised'

options, remainder = parser.parse_args()

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

print 'Starting ... '


# print N_processes
print folder_list


#using cluster dictionary 
if options.cluster_dictionary is not "None":  
       
    def plot_folder1(folder):
        if os.path.isdir(folder):
            print 'Plotting: ', folder
            folder_id = int(string.split(folder, '/')[-1])
            subprocess.call('source ~/.bashrc; python analyse_solutions_from_traces.py -p {0} -d {1} -c {2} -i {3}'.format(options.parfile,folder,options.cluster_dictionary,folder_id),
                            shell=True, executable='/bin/bash')
            subprocess.call('source ~/.bashrc; python taurex_mlehess.py -p {0} -d {1} -c {2} -i {3}'.format(options.parfile,folder,options.cluster_dictionary,folder_id),
                            shell=True, executable='/bin/bash')

    
#     plot_func1 = functools.partial(plot_folder1,params)
    
    
    pool = mp.Pool(processes=int(N_processes))           #setting number of cores on which to run
    pool_result = pool.map(plot_folder1,folder_list) #runnning the stuff
    pool.close()                                    #closing the pool 
        

#just indexing a folder and run things in there 
else:
    def plot_folder2(folder):
        if os.path.isdir(folder):
            print 'Plotting: ', folder
            if options.parfile is None:
                parfile = glob.glob(folder+'/*.par')[0]
            else:
                parfile = options.parfile
            subprocess.call('source ~/.bashrc; python analyse_solutions_from_traces.py -p '+parfile+' -d '+folder,shell=True,executable='/bin/bash')
#             os.system('python taurex_mlehess.py -p '+parfile+' -d '+folder)
            
    pool = mp.Pool(processes=int(N_processes))           #setting number of cores on which to run
    pool_result = pool.map(plot_folder2,folder_list) #runnning the stuff
    pool.close()                                    #closing the pool 
    # pool.join()                                     #closing the pool.

