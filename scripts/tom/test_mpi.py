
from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank==0:
    version= MPI.Get_version()
    sys.stdout.write("mpi version is {}\n".format(version))

    size = comm.Get_size()
    sys.stdout.write("size is {}\n".format(size))
