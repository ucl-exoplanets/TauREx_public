'''
    TauREx v2 - Development version - DO NOT DISTRIBUTE

    TauREx execution code

    Developers: Ingo Waldmann, Marco Rocchetto (University College London)

'''

# loading libraries
import  sys
import os
import argparse
import logging

# setting some parameters to global
global MPIimport, MPIverbose, MPImaster, multinest_import

#trying to initiate MPI parallelisation
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    MPIsize = MPI.COMM_WORLD.Get_size()
    if MPIsize > 1:
        MPIimport = True
        if MPIrank == 0:
            MPImaster = True
            MPIverbose= True
        else:
            MPImaster = False
            MPIverbose= False
    else:
        MPIverbose = True
        MPImaster  = True
        MPIimport  = False
except ImportError:
    MPIimport = False
    MPImaster = True
    MPIverbose= True


# checking for multinest library
try: 
    import pymultinest
    multinest_import = True
except:
    multinest_import = False

# loading classes
sys.path.append('./classes')
sys.path.append('./library')

from parameters import *
from library_constants import *

def node_architecture():
    """
    Check the presenze of GPU cards or cpu in the using node
    :return: boolean list for gpu and cpu accessible components
    """
    hard_list=[]
    # Check for an existing GPU first
    gpu_check = "lspci | grep -i 'vga\|3d\|2d'"
    gpu_list = os.popen(gpu_check).read()

    if 'Tesla K40c' in gpu_list:
        hard_list.append(True)
    else:
        hard_list.append(False)

    #Check for cpu cores
    cpu_check = "lscpu | egrep '^CPU\('"
    cpu_number = int(os.popen(cpu_check).read()[8:])
    # print "Available cpu cores: ", cpu_number

    hard_list.append([True]*cpu_number)
    return hard_list


processors = node_architecture()



#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('-p',
                    dest='param_filename',
                    default='Parfiles/default.par'
                   )

# cluster parameters
parser.add_argument('-c', '--cluster',
                  dest='cluster_dictionary',
                  default='None',
                  action='store'
                  )
parser.add_argument('-i', '--cluster_index',
                  dest='cluster_procid',
                  default='None',
                  type = str,
                  )

# misc  parameters
parser.add_argument('--no_multinest',
                  dest='no_multinest',
                  action='store_true',
                  default=False)

# plotting parameters
parser.add_argument('--plot',
                      action='store_true',
                      dest='plot',
                      default=False)
parser.add_argument('--plot_profiles',
                      action='store_true',
                      dest='plot_profiles',
                      default=False)
parser.add_argument('--plot_title',
                      dest='plot_title',
                      default=False)
parser.add_argument('--plot_prefix',
                      dest='plot_prefix',
                      default=False)
parser.add_argument('--plot_out_folder',
                      dest='plot_out_folder',
                      default=False)


# add command line interface to parameter file

params_tmp = parameters(log=False)
params_dict = params_tmp.params_to_dict() # get all param names
for param in params_dict:
    if type(params_dict[param]) == list:
        # manage lists
        parser.add_argument('--%s' % param,
                            action='append',
                            dest=param,
                            default=None,
                            type = type(params_dict[param][0])
                            )
    else:
        parser.add_argument('--%s' % param,
                            dest=param,
                            default=None,
                            type = type(params_dict[param])
                            )
options = parser.parse_args()

# Initialise parameters instance
params = parameters(options.param_filename, mode='retrieval')

params.gen_manual_waverange = False

# Override params object from command line input
for param in params_dict:
    if getattr(options, param) != None:
        value = getattr(options, param)
        if param == 'planet_mass':
            value *= MJUP
        if param == 'planet_radius':
            value *= RJUP
        if param == 'star_radius':
            value *= RSOL
        if param == 'atm_mu':
            value *= AMU
        setattr(params, param, value)


# MPI related message
if MPIimport:
    logging.info('MPI enabled. Running on %i cores' % MPIsize)
else:
    logging.warning('MPI disabled')

if multinest_import:
    logging.info('Multinest library correctly loaded.')
else:
    logging.warning('Multinest library is not loaded. Multinest functionality will be disabled')

# modifying parameters object if running in cluster mode (see cluster class for docu)
if options.cluster_dictionary is not 'None':
    from cluster import cluster
    c = cluster()
    c.read_dict(dict_name=options.cluster_dictionary)
    params = c.modify_params(params,options.cluster_procid)



if params.gen_type == 'transmission':
    from taurex_transmission import run
elif params.gen_type == 'emission':
    from taurex_emission import run
else:
    logging.error('Forward model selected is ambiguous')
    logging.info('Check \'type\' and \'fit_emission\', \'fit_transmission\' parameters')
    logging.info('PS: you suck at this... ')
    exit()

if MPIimport:
    MPI.COMM_WORLD.Barrier() # wait for everybody to synchronize here

#running Tau-REx
if MPIimport:
    if processors[0]:
        if MPI.COMM_WORLD.Get_rank() == 0:
            outputob = run(params, options, gpu=True)
        else:
            outputob = run(params, options, gpu=False)
    else:
        outputob = run(params, options, gpu=False)

# plotting
if  options.plot:
    sys.path.append('./tools')
    from taurex_plots import taurex_plots
    for val in outputob.dbfilename:
        logging.info('Initialising plotting')
        dbfilename = outputob.dbfilename[val]
        plot = taurex_plots(pickle_fname=dbfilename, title=options.plot_title, prefix=options.plot_prefix,
                            out_folder=options.plot_out_folder)
        logging.info('Plot posterior distributions')
        plot.plot_posteriors()
        logging.info('Plot fitted spectrum')
        plot.plot_fitted_spectrum()
        if params.gen_type == 'emission' or params.gen_ace or options.plot_profiles:
            logging.info('Plot mixing ratio profiles and temperature pressure profile')
            plot.plot_fitted_xprofiles()
            plot.plot_fitted_tp()