#!/bin/sh
#PBS -q celltypes
#PBS -N mpi_test
#PBS -r n
#PBS -l walltime=00:01:00

cd $PBS_O_WORKDIR

source activate aaopt3
echo "Conda MPI"
echo "PBS_np: ($PBS_NP)"
echo $(which mpiexec)
# mpiexec python /home/tom.chartrand/work/ateam-tools/scripts/tom/test_mpi.py -outfile-pattern engine1_%r-%g-%h -output-filename test_out &
mpiexec -outfile-pattern engine1_%r python -c "print('test')" &
sleep 3
module add mpi/mpich-3.2-x86_64
echo "MPI  3.2"
echo "PBS_np: ($PBS_NP)"
echo $(which mpiexec)
# mpiexec python /home/tom.chartrand/work/ateam-tools/scripts/tom/test_mpi.py -outfile-pattern engine2_%r-%g-%h -output-filename test_out &
mpiexec -outfile-pattern engine2_%r python -c "print('test')" &
sleep 3