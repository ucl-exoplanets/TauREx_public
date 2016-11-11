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
        
    
    def generate_script(self, ID_number,DICT=None,DICTFILENAME = None):
        '''
        function generating the submit script. This can be improved (I guess) 
        but works for now.
        '''
  
        if DICTFILENAME is None:
            DICTFILE = self.DICTNAME 
        if DICT is None:
            DICT = self.DICT
        
       
        #getting parameters
        NODES = DICT[str(ID_number)]['GENERAL']['NODES']
        CPUS  = DICT[str(ID_number)]['GENERAL']['CPUS']
        WTIME = DICT[str(ID_number)]['GENERAL']['WALLTIME']
        MEM   = DICT[str(ID_number)]['GENERAL']['MEMORY']
        PFILE = DICT[str(ID_number)]['GENERAL']['PARFILE']
        DISKMEM = DICT[str(ID_number)]['GENERAL']['DISKMEM']
        CLUSTER = DICT[str(ID_number)]['GENERAL']['CLUSTER']
        USERID  = DICT[str(ID_number)]['GENERAL']['USERNAME']
        TCPUS = int(NODES)*int(CPUS) #total number of cores used
        
        if CLUSTER is 'cobweb':
            SCRATCHDIR = DICT[str(ID_number)]['GENERAL']['SCRATCH_DIR']
        OUTDIR  = DICT[str(ID_number)]['GENERAL']['OUTPUT_DIR']
        DATADIR = DICT[str(ID_number)]['GENERAL']['DATA_DIR']
        
#         print 'submit_{}.sh'.format(ID_number)
#         exit()
        scriptname = 'submit_{}.sh'.format(ID_number)        
        if CLUSTER is 'cobweb':
            self.__generate_cobweb(ID_number, scriptname, TCPUS, PFILE, DICTFILE, NODES, CPUS, WTIME, MEM, SCRATCHDIR, DATADIR, OUTDIR)
        elif CLUSTER is 'legion':
            self.__generate_legion(ID_number,scriptname,USERID,TCPUS,PFILE,DICTFILE,WTIME,MEM,DISKMEM,DATADIR,OUTDIR)
    
        return scriptname
    
    def __generate_cobweb(self,ID_number,scriptname,TCPUS,PFILE,DICTFILE,NODES,CPUS,WTIME,MEM,SCRATCHDIR,DATADIR,OUTDIR):
            with open(scriptname,'wb') as h:
                #general setup stuff for the COBWEB cluster
                
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
                h.write('rm -rf {0}/{1}'.format(SCRATCHDIR,ID_number)+'\n')                     #removes previous stuff if there
                h.write('mkdir -p {0}/{1}'.format(SCRATCHDIR,ID_number)+'\n')                   #make tmp folder on scratch
                h.write('cp -rf {0}/* {1}/{2}'.format(DATADIR,SCRATCHDIR,ID_number)+'\n')       #copying everything over to scratch
                h.write('cd {0}/{1}'.format(SCRATCHDIR,ID_number)+'\n')                         #we are going places
                h.write('cat $PBS_NODEFILE > nodes'+'\n')              #creating hostfile 
                
                #run command standard version
                h.write('mpirun -np {0} -hostfile nodes  /share/apps/anaconda/2.2.0/bin/python taurex.py -p {1} -c {2} -i {3}'.format(TCPUS,PFILE,DICTFILE,ID_number)+' --plot \n') 
                
                #run command version for multiple processes per node (i.e. ppn <24)
    #             h.write('OMPI_MCA_btl=^openib OMPI_MCA_mtl=^psm mpirun -np {0} -hostfile nodes  /share/apps/anaconda/2.2.0/bin/python taurex.py -p {1} -c {2} -i {3}'.format(TCPUS,PFILE,DICFILE,ID_number)+'\n') #run command
    
                #version without hostfile. works better if run on 3+ nodes 
    #             h.write('mpirun -np {0} /share/apps/anaconda/2.2.0/bin/python taurex.py -p {1} -c {2} -i {3}'.format(TCPUS,PFILE,DICFILE,ID_number)+'\n') #run command
    
                h.write('mkdir -p {0}'.format(OUTDIR)+'\n')                       #make output dir
                h.write('cp -rf {0}/{1}/Output/* {2}'.format(SCRATCHDIR,ID_number,OUTDIR)+'\n')    #copy all the stuff back
                h.write('rm -rf {0}/{1}'.format(SCRATCHDIR,ID_number)+' \n')                    #remove scratch dir
    
    
    
    def __generate_legion(self,ID_number,scriptname,USERID,TCPUS,PFILE,DICTFILE,WTIME,MEM,DISKMEM,DATADIR,OUTDIR):
            #general setup stuff for the LEGION cluster 
        with open(scriptname,'wb') as h:
            h.write('#! /bin/bash -l'+'\n')                      #sets script shell
            h.write('#$ -S /bin/bash'+'\n')                      #sets runtime shell
            h.write('#$ -N taurex_run_{}'.format(ID_number)+'\n')#sets runtime name 
            h.write('#$ -V'+'\n')                                #imports login-node environment & path variables
            
            #hardware setup
            h.write('#$ -l h_rt={}'.format(WTIME)+'\n')          #total runtime
            h.write('#$ -l mem={}G'.format(MEM)+'\n')           #total memory used 
            h.write('#$ -l tmpfs={}G'.format(DISKMEM)+' \n')    #total amount of disk space needed (default 10GB)   
            h.write('#$ -pe mpi {}'.format(TCPUS)+'\n')          #total number of cpus 
            
            #setting working directory 
            h.write('#$ -wd '+DATADIR+'\n') #setting up path to SCRATCH space, data needs to be in scratch not home 
            
            #making sure environment variables are set
#             h.write('module unload compilers/intel/2015/update2 \n')
#             h.write('module unload mpi/intel/2015/update3/intel \n')
#             h.write('module load compilers/gnu/4.9.2 \n')
#             h.write('module load mpi/openmpi/1.10.1/gnu-4.9.2 \n')
#             h.write('module load mpi4py/2.0.0/python2 \n')
            h.write('export MPLBACKEND="pdf" \n')
            
            #setting up output dir on scratch
            outdirpath = OUTDIR#+'/{0}'.format(ID_number)
            h.write('mkdir -p '+outdirpath+'\n')  #setting up output directory
        
            #run main command 
            h.write('gerun python taurex.py -p {0} -c {1} -i {2}'.format(PFILE,DICTFILE,ID_number)+' --plot \n') 
    

         
            
        
            
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
            print('Error: dictionary file found in base directory.')
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
        
        
        
        
    
    
    
    
    