# loading libraries
import  sys
import os


#trying to initiate MPI parallelisation
try:
    from mpi4py import MPI
    MPIrank = MPI.COMM_WORLD.Get_rank()
    MPIsize = MPI.COMM_WORLD.Get_size()
    MPIimport = True
except ImportError:
    MPIimport = False

# checking for multinest library
try: 
    import pymultinest
    multinest_import = True
except:
    multinest_import = False


if MPIimport: 
    print 'MPI imported'
else:
    print 'no MPI'
    
if multinest_import:
    print 'imported pymultinest, all is well'
else: 
    print 'not imported pymultinest'