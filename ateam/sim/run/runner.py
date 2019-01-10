"""Functions for running bionet simulations, including with MPI and submitting batch jobs."""

from bmtk.simulator import bionet
from ateam.sim import cell_functions 
import subprocess
import os

_pycommand = "from ateam.sim.run import runner; runner.run_bionet(\\\"{config}\\\")"

def bionet_mpi_command(config, ncores=1):
    command = _pycommand.format(config=config)
    mpicommand = "mpirun -np {ncores} nrniv -mpi -python -c \"{pycommand}\"".format(ncores=ncores, pycommand=command)
    return mpicommand

def run_bionet_mpi(config, ncores=1):
    mpicommand = bionet_mpi_command(config, ncores)
    sp = subprocess.Popen([mpicommand], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    # error first
    print out[1] 
    print out[0]

def run_bionet(config):
    # os.chdir(os.path.dirname(config))
    conf = bionet.Config.from_json(config, validate=True)
    conf.build_env()
    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    return sim
    # bionet.nrn.quit_execution()

def run_bionet_hpc(config):
    pass