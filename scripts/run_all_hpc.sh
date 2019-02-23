#!/bin/sh


for i in $(find $1 -type f -path '*/psp/config.json'); 
do
    run_hpc "$i" -t 00:05:00
done
for i in $(find $1 -type f -path '*/rate/config.json'); 
do
    run_hpc "$i" -t 00:15:00 -n 2
done