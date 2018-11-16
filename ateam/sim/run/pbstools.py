__version__ = '0.1.0'

import argparse
import copy
import subprocess as sp
import logging
import os



logger = logging.getLogger(__name__)

default_settings = {
                    'jobname':'Default_Job_Name',
                    'email':None,
                    'email_options':None,
                    'queue':'celltypes',
                    'mem':'4g',
                    'vmem':None,
                    'walltime':'01:00:00',
                    'ncpus':None,
                    'ppn':None,
                    'nodes':None,
                    'jobdir':None,
                    'outfile':'$PBS_JOBID.out',
                    'errfile':'$PBS_JOBID.err',
                    'environment':{},
                    'array':None,
                    'priority':'low',
                    'rerunable':False
                    }

class PBSJob(object):

    qsub_command = 'qsub'

    def __init__(self, command, **kwargs):

        settings_dict = copy.copy(default_settings)
        settings_dict.update(kwargs)
        for key, val in settings_dict.iteritems():
            setattr(self, key, val)

        self.command = command

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

        script_lines = []

        script_lines.append('#!/bin/bash\n')

        if not self.array is None:
            script_lines.append('#PBS -t %s-%s\n' % (self.array[0], self.array[1]))

        if not self.priority is None:
            try:
                assert self.priority in ['high', 'med', 'low']
            except AssertionError:
                raise ValueError('Priority "%s" is not "high/med/low"')

            script_lines.append('#PBS -W x=QOS:%s\n' % (self.priority,))

        script_lines.append('#PBS -q %s\n' % (self.queue))
        script_lines.append('#PBS -N %s\n' % (self.jobname))

        if self.email != None:
            script_lines.append('#PBS -M %s\n' % (self.email))

        if self.email_options != None:
            script_lines.append('#PBS -m %s\n' % (self.email_options))

        if self.rerunable == False:
            script_lines.append('#PBS -r n\n')

        if not self.procs is None:
            script_lines.append('#PBS -l procs=%d\n' % self.procs)

        cpu_ppn_node_list = []

        if not self.ncpus is None:
            cpu_ppn_node_list.append('ncpus=%d' % self.ncpus)

        if not self.nodes is None:
            cpu_ppn_node_list.append('nodes=%s' % self.nodes)

        if not self.ppn is None:
            cpu_ppn_node_list.append('ppn=%d' % self.ppn)

        if len(cpu_ppn_node_list) > 0:
            tmp_str = '#PBS -l %s\n' % ':'.join(cpu_ppn_node_list)
            script_lines.append(tmp_str)

        vmem_walltime = []
        if self.mem != None:
            vmem_walltime.append('mem=%s' % (self.mem))

        if self.vmem != None:
            vmem_walltime.append('vmem=%s' % (self.vmem))

        if self.walltime != None:
            vmem_walltime.append('walltime=%s' % (self.walltime))

        if len(vmem_walltime) > 0:
            tmp_str = '#PBS -l %s\n' % (','.join(vmem_walltime))
            script_lines.append(tmp_str)

        if self.jobdir != None:
            script_lines.append('#PBS -d %s\n' % (os.path.expanduser(self.jobdir)))

        script_lines.append('#PBS -o %s\n' % (os.path.expanduser(self.outfile)))
        script_lines.append('#PBS -e %s\n' % (os.path.expanduser(self.errfile)))

        env_list = []
        for variable, value in self.environment.items():
            env_list.append('%s=%s' % (variable, value))

        if len(env_list) > 0:
            script_lines.append('#PBS -v %s \n' % (','.join(env_list)))

        if self.conda_env != None:
            script_lines.append('source activate %s\n' % self.conda_env )

        script_lines.append('%s\n' % (self.command))

        script_string = ''.join(script_lines)
        logger.info(script_string)

        for line in script_lines:
            if verbose: print line,
            if not dryrun: sp_input.write(line)

        if not dryrun:
            sp_input.close()
            result = sp_output.read()
            if verbose: print result
        else:
            result = None
        return result

class PythonJob(PBSJob):

    def __init__(self, script=None, command=None, python_executable='python', python_args='', python_path = None, **kwargs):

        self.python_executable = python_executable
        # self.conda_env = conda_env
        # assert os.path.exists(script)
        assert (script or command) and not (script and command)
        self.pycommand = "-c \"{command}\"".format(command=command) if command else script
        self.python_args = python_args

        command = '%s %s %s' % (self.python_executable, self.pycommand, self.python_args)

        # if not conda_env is None:
        #     command = 'source activate %s\n' % self.conda_env + command

        if not python_path is None:
            command = ("export PYTHONPATH=%s\n" % python_path) + command

        super(PythonJob, self).__init__(command, **kwargs)

pbs_parser = argparse.ArgumentParser()
for key, val in default_settings.iteritems():
    pbs_parser.add_argument('--%s' % key, default=str(val), type=str)



if __name__ == "__main__":  # pragma: no cover
    pass
