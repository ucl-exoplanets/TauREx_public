#! /usr/bin/python -W ignore

###########################################################
# cluster class: 
#
# This class both generates qsub submit.sh scripts automatically 
# submits these scripts and dynamically modifies the taurex params 
# object on runtime
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, June 2014      
#       
###########################################################


#loading libraries     
import sys, os, optparse, time,json
import numpy as np #nummerical array library 


#loading classes
# sys.path.append('./classes')
# sys.path.append('./library')



class cluster(object):
    '''
    The cluster object does two things: 1) it generates and manages submit scripts;
    2) it modifies the params object on runtime. 
    Before running, cluster object requires an ascii file with a parsable python dictionary following 
    this format: 
    
    input_dict[ID_number][parameter1][value]     #for parameter object parameters to be modified on runtime
    input_dict[ID_number][parameter2][value]
    
    input_dict[ID_number][GENERAL][NODES]       #general run parameters
                                  [CPUS]
                                  [WALLTIME]    #hh:mm:ss
                                  [MEMORY]      #in GB
                                  [PARFILE]
                                  [OUTPUT_DIR]
                                  
    1) The cluster class will generate a qsub submit_$ID_number.sh script for each ID_number and run taurex with 
    
    #PBS -l nodes=$NODES:ppn=$CPUS
    mpirun -np $CPUs python taurex.py -p $PARFILE -c $ID_number
    cp -rf /scratch/run/$ID_number/* $OUTPUT_DIR_PATH
    
    2) On runtime taurex will read the same input_dict and modify params in parameter object for given ID_number
    '''

    def __init__(self,options=None):
        '''
        options: command line parsed options
                  -i input_dict name
                  -d data directory (/share/data/...)
        '''
        if options is not None:
            self.DICTNAME = options.dictname
            self.DATADIR  = options.datadir
        
    
    def generate_script(self, ID_number,DICT=None,DICTFILENAME = None, DATADIR=None):
        '''
        function generating the submit script. This can be improved (I guess) 
        but works for now.
        '''
        SCRATCHDIR = '/scratch/ingo/run/{}'.format(ID_number)
        if DATADIR is None:
            DATADIR    = self.DATADIR
        if DATADIR[-1] is '/':
            DATADIR = DATADIR[:-1]
        if DICTFILENAME is None:
            DICFILE = self.DICTNAME 
        if DICT is None:
            DICT = self.DICT
        
       
        #getting parameters
        NODES = DICT[str(ID_number)]['GENERAL']['NODES']
        CPUS  = DICT[str(ID_number)]['GENERAL']['CPUS']
        WTIME = DICT[str(ID_number)]['GENERAL']['WALLTIME']
        MEM   = DICT[str(ID_number)]['GENERAL']['MEMORY']
        PFILE = DICT[str(ID_number)]['GENERAL']['PARFILE']
        OUTDIR= DICT[str(ID_number)]['GENERAL']['OUTPUT_DIR']
        TCPUS = int(NODES)*int(CPUS) #total number of cores used
        
#         print 'submit_{}.sh'.format(ID_number)
#         exit()
        
        scriptname = 'submit_{}.sh'.format(ID_number)
        
        with open(scriptname,'wb') as h:
            #general setup stuff
            h.write('#! /bin/bash'+'\n')                           #sets script shell
            h.write('#PBS -S /bin/bash'+'\n')                      #sets runtime shell
            h.write('#PBS -N taurex_run_{}'.format(ID_number)+'\n')#sets runtime name 
            h.write('#PBS -q compute'+'\n')                        #sets queue: compute/test
            h.write('#PBS -j oe'+'\n')                             #merges STDOUT with SDTERR
            #h.write('#PBS -o ./run_log_%s.txt',(ID_number))       #writes STDOUT to specific file
            h.write('#PBS -V'+'\n')                                #imports login-node environment & path variables
            
            #hardware setup
            h.write('#PBS -l nodes={}:ppn={}'.format(NODES,CPUS)+'\n')   #number of nodes and cores per node
            h.write('#PBS -l walltime={}'.format(WTIME)+'\n')            #total runtime
            h.write('#PBS -l mem={}gb'.format(MEM)+'\n')                 #total memory used (max 256GB per node)
            
            #copy and run
            h.write('mkdir -p '+SCRATCHDIR+'\n')                   #make tmp folder on scratch
            h.write('cp -rf '+DATADIR+'/* '+SCRATCHDIR+'\n')       #copying everything over to scratch
            h.write('cd '+SCRATCHDIR+'\n')                         #we are going places
            h.write('cat $PBS_NODEFILE > nodes'+'\n')              #creating hostfile 
            h.write('mpirun -np {0} -machinefile nodes  /share/apps/anaconda/2.2.0/bin/python taurex.py -p {1} -c {2} -i {3}'.format(TCPUS,PFILE,DICFILE,ID_number)+'\n') #run command
            h.write('mkdir -p '+OUTDIR+'\n')                       #make output dir
            h.write('cp -rf '+SCRATCHDIR+'/* '+OUTDIR+'/'+'\n')    #copy all the stuff back
            
        return scriptname
            
    def write_dict(self,DICT,DICTNAME='run_dict.dat'):
        '''
        function that takes externally generated dictionary and dumps it to ascii.
        Using this function avoids compatibility issues later on.
        '''
        dicdump = json.dumps(DICT)
        with open(DICTNAME,'wb') as outfile:
            outfile.write(dicdump+'\n')
    
       
    def read_dict(self,dict_name=None):
        '''
        function reading in an parsing ascii dictionary. It sets the self.DICT variable. 
        '''
        if dict_name is None:
            dict_name = self.DICTNAME
            
        try:
            with open(dict_name,'r') as infile:
                DICT_in = infile.readlines()
        except IOError:
            print 'Error: dictionary file found in base directory.'
            exit()
                
        #parsing dictionary
        self.DICT = eval(DICT_in[0],{'false': False, 'true': True, 'null': None})
        self.IDs  = self.DICT.keys()
    

    
    def modify_params(self,params,ID_number):
        '''
        function replacing parameters in params with dictionary values
        '''
        keys = self.DICT[str(ID_number)].keys()
#         N_keys = len(keys)
    
        for par in keys:
            setattr(params,par, self.DICT[str(ID_number)][par])
            
        return params
        
        
        
        
    
    
    
    
    