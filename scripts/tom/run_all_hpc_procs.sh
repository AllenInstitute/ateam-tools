#!/bin/sh

DIR=${1-'.'}
for i in $(find $1 -type f -path "$DIR/*/psp/config.json"); 
do
    run_hpc_bmtk "$i" -t 00:30:00 -p 64 --jobname hcem_singlecell_sim --conda bmtk
done