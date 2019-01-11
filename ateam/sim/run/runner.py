"""Functions for running bionet simulations, including with MPI and submitting batch jobs."""

import subprocess
import os

# TODO: could put cd to config dir in subproc commands?
_pycommand = r'from ateam.sim.run import run_bionet; run_bionet.run(\"{config}\")'

def bionet_mpi_command(config, ncores=1):
    command = _pycommand.format(config=config)
    mpicommand = 'mpirun -np {ncores} nrniv -mpi -python -c "{pycommand}"'.format(ncores=ncores, pycommand=command)
    return mpicommand

def run_bionet_mpi(config, ncores=1):
    mpicommand = bionet_mpi_command(config, ncores)
    sp = subprocess.Popen([mpicommand], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    # error first
    print out[1] 
    print out[0]

def run_bionet(config):
    command = 'python -c "{cmd}"'.format(cmd=_pycommand).format(config=config)
    sp = subprocess.Popen([command], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    # error first
    print out[1] 
    print out[0]