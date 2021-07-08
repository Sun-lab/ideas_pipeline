#!/bin/bash

# sbatch -c 8 --mem 32G step1b_MAST_l2.sh
# srun -N 1 -n 1 -c 4 --mem=32G --pty bash
# sstat -a -j 31623384,31515922 -o jobid,averss,maxrss,avevmsize,maxvmsize
# ml R/4.1.0-foss-2020b
for type in $( cat cell_types_mast_l2.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_MAST.R step1b_Rout/step1b_MAST_${type}.Rout
done
