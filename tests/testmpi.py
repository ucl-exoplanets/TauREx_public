#!/usr/bin/env python

from mpi4py import MPI
import thread
import numpy as np
import os

comm = MPI.COMM_WORLD
 
#print 'groupsize ', MPI.Group_size(comm.Get_group())
 
print "Hello! I'm rank %d from %d in %s" % (comm.rank, comm.size,comm.name )
 
comm.Barrier()   # wait for everybody to synchronize _here_

test = np.ones((1,100))

if comm.rank == 0:
    print 'Writing test file'
    np.savetxt('testfile.txt',test)
    print 'write done'
    print 'Writing test folder'
    os.mkdir('testdir2')
    print 'testdir2 done'


#comm=MPI.COMM_WORLD
#rank=comm.Get_rank()
#size=comm.Get_size()
#pares=[0,1,2,3]
#pares2 = [4,5,6,7]
#
#group=comm.Get_group()
#newgroup=MPI.Group.Incl(group,pares)
#newcomm=comm.Create(newgroup)
#if rank < 4:
#    print newcomm.Get_rank(), newcomm.Get_size()
#
## newcomm.Barrier()
#
#newgroup2 = MPI.Group.Incl(group,pares2)
#newcomm2 = comm.Create(newgroup2)
#if rank > 4:
#    print newcomm.Get_rank(), newcomm.Get_size()
## print newcomm2.Get_rank(),newcomm2.Get_size()
#
#
#comm.Barrier()

# newgroup.Free()
# group.Free()
# newcomm.Free()
