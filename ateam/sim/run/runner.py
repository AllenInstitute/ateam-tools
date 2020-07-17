"""Functions for running bionet simulations, including with MPI and submitting batch jobs."""

import subprocess
import os
import os.path

_pycommand = r'from ateam.sim.run import run_bionet; run_bionet.run(\"{config}\")'

def bionet_mpi_command(config, ncores=1):
    command = _pycommand.format(config=config)
    sim_dir = os.path.dirname(config)
    mpicommand = 'cd {dir}; mpirun -np {ncores} nrniv -mpi -python -c "{pycommand}"'.format(dir=sim_dir, ncores=ncores, pycommand=command)
    return mpicommand

def run_bionet_mpi(config, ncores=1):
    mpicommand = bionet_mpi_command(config, ncores)
    sp = subprocess.Popen([mpicommand], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    # error first
    print(out[1], out[0])

def run_bionet(config):
    sim_dir = os.path.dirname(config)
    command = 'cd {dir}; python -c "{cmd}"'.format(dir=sim_dir, cmd=_pycommand).format(config=config)
    sp = subprocess.Popen([command], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    # error first
    print(out[1], out[0])