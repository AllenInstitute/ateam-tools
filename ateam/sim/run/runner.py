
from bmtk.simulator import bionet

pycommand = "from ateam.sim.run import runner; runner.run_bionet(\\\"{config}\\\")"

def run_bionet_mpi(config, ncores=1):
    import subprocess
    # config = os.path.basename(config)
    command = pycommand.format(config=config)
    mpicommand = "mpirun -np {ncores} nrniv -mpi -python -c \"{pycommand}\"".format(ncores=ncores)
    
    sp = subprocess.Popen([mpicommand], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    out = sp.communicate()
    print out[0]; print out[1]


def run_bionet(config):
    conf = bionet.Config.from_json(config, validate=True)
    conf.build_env()
    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    # bionet.nrn.quit_execution()

def run_bionet_hpc(config):
    pass