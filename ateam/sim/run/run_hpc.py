'''Script for job submission on HPC
Generates a PBS job submission script from the given options.
Activates a Conda environment named 'bmtk' when starting the job
'''

from ateam.sim.run import pbstools
from ateam.sim.run import runner
import os.path
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("config")
parser.add_argument("--conda", default="bmtk_ateam", help="Conda environment name")
parser.add_argument("--time","-t", default="00:10:00")
parser.add_argument("--nodes", "-n", default=1, type=int)
parser.add_argument("--mem", "-m", default="4g", help="Memory usage (only used for single node jobs?)")
parser.add_argument("--procs","-p", type=int, help="Number of processors. Overrides nodes and ppn.")
parser.add_argument("--queue","-q", default="celltypes")
parser.add_argument("--ppn", type=int)
parser.add_argument("--norerun", action='store_true', help="Prevent overwriting existing sim output")
parser.add_argument("--export_vars", action='store_true')

def main(args_list=None):
    args = parser.parse_args(args_list)
    if args.procs:
        args.nodes = None
        args.ppn = None
        n = args.procs
    else:
        ppn = args.ppn or 24
        n = ppn*args.nodes


    options = {
        'queue':args.queue,
        'jobname':'bmtk_test',
        'walltime':args.time,
        'nodes':args.nodes,
        'ppn':args.ppn,
        'mem':args.mem,
        'procs':args.procs,
        'priority':None
    }

    folder = os.path.dirname(args.config)
    if args.norerun and os.path.isfile(os.path.join(folder, "output", "spikes.h5")):
        print("Output already exists, aborting.")
        return
    command = runner.bionet_mpi_command(config=args.config, ncores=n)
    job = pbstools.PBSJob(command=command, jobdir=folder, conda_env=args.conda, **options)
    job.run()
    