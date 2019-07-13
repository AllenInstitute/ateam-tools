#!/bin/sh


for i in $(find $1 -type f -path '*/psp/config.json'); 
do
    run_hpc "$i" -t 00:10:00 --norerun --jobname hcem_singlecell_sim
done
for i in $(find $1 -type f -path '*/rate/config.json'); 
do
    run_hpc "$i" -t 00:20:00 -n 4 --norerun --jobname hcem_singlecell_sim
done