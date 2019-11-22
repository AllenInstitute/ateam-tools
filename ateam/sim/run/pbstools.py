import argparse
import copy
import subprocess as sp
import os

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("--conda", default="bmtk_ateam", help="Conda environment name")
parser.add_argument("--walltime","-t", default="00:10:00")
parser.add_argument("--nodes", "-n", default=1, type=int)
parser.add_argument("--mem", "-m", default="4g", help="Memory usage (only used for single node jobs?)")
parser.add_argument("--procs","-p", type=int, help="Number of processors. Overrides nodes and ppn.")
parser.add_argument("--queue","-q", default="celltypes")
parser.add_argument("--jobname", default="ateam_test_job")
parser.add_argument("--jobdir", default=os.getcwd())
parser.add_argument("--ppn", type=int)
parser.add_argument("--ncpus", type=int)
parser.add_argument("--vmem", default=None)
parser.add_argument("--email", default=None)
parser.add_argument("--email_options", default='a')
parser.add_argument("--array", default=None)
parser.add_argument("--environment", default={})
parser.add_argument("--outfile", default='$PBS_JOBID.out')
parser.add_argument("--errfile", default='$PBS_JOBID.err')
parser.add_argument("--rerun", action='store_true')
parser.add_argument("--dryrun", action='store_true')
# subparsers = parser.add_subparsers()


class PBSJob(object):

    qsub_command = 'qsub'

    def __init__(self, command, args):
        if args.procs:
            args.nodes = None
            args.ppn = None
            self.n = args.procs
        else:
            ppn = args.ppn or 24
            self.n = ppn*args.nodes
        self.command = command
        self.script_lines = self.generate_script(args)

    def generate_script(self, args):
        script_lines = []
        script_lines.append('#!/bin/bash\n')
        script_lines.append('#PBS -q %s\n' % (args.queue))
        script_lines.append('#PBS -N %s\n' % (args.jobname))

        if args.array is not None:
            script_lines.append('#PBS -t %s-%s\n' % (args.array[0], args.array[1]))
        if args.email is not None:
            script_lines.append('#PBS -M %s\n' % (args.email))
        script_lines.append('#PBS -m %s\n' % (args.email_options))
        script_lines.append('#PBS -r {}\n'.format('y' if args.rerun else 'n'))

        if args.procs is not None:
            script_lines.append('#PBS -l procs={}\n'.format(args.procs))

        cpu_ppn_node_list = []
        # if not args.ncpus is None:
        #     cpu_ppn_node_list.append('ncpus=%d' % args.ncpus)
        if args.nodes is not None:
            cpu_ppn_node_list.append('nodes=%s' % args.nodes)
        if args.ppn is not None:
            cpu_ppn_node_list.append('ppn=%d' % args.ppn)
        elif not args.procs:
            script_lines.append('#PBS -n \n')

        if len(cpu_ppn_node_list) > 0:
            tmp_str = '#PBS -l %s\n' % ':'.join(cpu_ppn_node_list)
            script_lines.append(tmp_str)

        vmem_walltime = []
        if args.mem != None:
            vmem_walltime.append('mem=%s' % (args.mem))
        if args.vmem != None:
            vmem_walltime.append('vmem=%s' % (args.vmem))
        if args.walltime != None:
            vmem_walltime.append('walltime=%s' % (args.walltime))

        if len(vmem_walltime) > 0:
            tmp_str = '#PBS -l %s\n' % (','.join(vmem_walltime))
            script_lines.append(tmp_str)

        if args.jobdir != None:
            script_lines.append('#PBS -d %s\n' % (os.path.expanduser(args.jobdir)))

        script_lines.append('#PBS -o %s\n' % (os.path.expanduser(args.outfile)))
        script_lines.append('#PBS -e %s\n' % (os.path.expanduser(args.errfile)))
        script_lines.append('#PBS -j oe\n')

        env_list = []
        for variable, value in args.environment.items():
            env_list.append('%s=%s' % (variable, value))

        if len(env_list) > 0:
            script_lines.append('#PBS -v %s \n' % (','.join(env_list)))

        if args.conda != None:
            script_lines.append('source activate %s\n' % args.conda )

        script_lines.append('%s\n' % (self.command))
        return script_lines

    def run(self, verbose=True, dryrun=False):

        if not dryrun:

            sub_process = sp.Popen(PBSJob.qsub_command,
                                   shell=True,
                                   stdin=sp.PIPE,
                                   stdout=sp.PIPE,
                                   close_fds=True)
            sp_output = sub_process.stdout
            sp_input = sub_process.stdin

            if sp_input == None:
                raise Exception('could not start job')

        for line in self.script_lines:
            if verbose: print line,
            if not dryrun: sp_input.write(line)

        if not dryrun:
            sp_input.close()
            result = sp_output.read()
            if verbose: print result
        else:
            result = None
        return result

class BmtkJob(PBSJob):
    def __init__(self, args):
        from ateam.sim.run import runner
        super(BmtkJob, self).__init__(None, args)
        args.jobdir = args.jobdir or os.path.dirname(args.config)
        self.command = runner.bionet_mpi_command(config=args.config, ncores=self.n)
        self.script_lines = self.generate_script(args)

def run_hpc(args_list=None):
    # parser = subparsers.add_parser()
    parser.add_argument("command")
    args = parser.parse_args(args_list)
    job = PBSJob(args.command, args)
    job.run(dryrun=args.dryrun)

def run_hpc_bmtk(args_list=None):
    parser.add_argument("config")
    parser.add_argument("--overwrite", action='store_true', help="Overwrite existing sim output")
    parser.set_defaults(jobname="bmtk_sim_test_job")
    parser.set_defaults(jobdir=None)
    args = parser.parse_args(args_list)
    # TODO: actually check output from config
    folder = os.path.dirname(args.config)
    if not args.dryrun and not args.overwrite and os.path.isfile(os.path.join(folder, "output", "spikes.h5")):
        print("Output already exists, aborting.")
        return
    job = BmtkJob(args)
    job.run(dryrun=args.dryrun)
