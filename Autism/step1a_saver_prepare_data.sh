#!/bin/sh

# sbatch -c 12 --mem 32G step1a_saver_prepare_data.sh
# sstat -j 31515346 -o jobid,averss,maxrss,avevmsize,maxvmsize

# ml R/4.0.3-foss-2020b

R CMD BATCH --no-save --no-restore step1a_saver_prepare_data.R step1a_saver_prepare_data.Rout

