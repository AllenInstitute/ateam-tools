#!/bin/sh


for i in $(find $1 -type f -path '*/psp/config.json'); 
do
    run_hpc "$i" -t 00:10:00 --jobname hcem_singlecell_sim --conda bmtk
done
for i in $(find $1 -type f -path '*/rate/config.json'); 
do
    run_hpc "$i" -t 00:20:00 -n 4 --jobname hcem_singlecell_sim --conda bmtk
done