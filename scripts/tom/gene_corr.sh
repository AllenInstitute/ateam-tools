#!/bin/sh
#PBS -q celltypes
#PBS -N genes_ephys
#PBS -r n
#PBS -l walltime=03:00:00
#PBS -l nodes=1:ppn=24

source activate ipfx3_home
python /home/tom.chartrand/work/ateam-tools/scripts/tom/gene_corr_depth_frem.py $PBS_ARRAYID
