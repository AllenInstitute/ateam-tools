#!/bin/sh


for i in $(find $1 -type f -path '*/psp/config.json'); 
do
    run_hpc_bmtk "$i" -t 00:10:00 -p 32 --jobname hcem_singlecell_sim --conda bmtk
done
# for i in $(find $1 -type f -path '*/rate/config.json'); 
# do
#     run_hpc_bmtk "$i" -t 00:10:00 -p 100 ---jobname hcem_singlecell_sim
# done