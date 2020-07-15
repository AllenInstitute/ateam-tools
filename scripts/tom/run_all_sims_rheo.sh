#!/bin/sh
source activate bmtk
DIR=${1-'.'}
for i in $(find $1 -type f -regex "$DIR/[0-9]+_kt/config.json"); 
do
    run_hpc_bmtk "$i" -t 00:40:00 -p 20 --jobname hcem_singlecell_sim --conda bmtk
done
