#!/bin/bash
source /app/lmod/lmod/init/profile

#ml R
export OMP_NUM_THREADS=1

sbatch -t 02-10:00:00 -c 32 \
--mem=128000 -o 1a_saver.log \
-p campus-new --job-name=1a_saver \
/fh/scratch/delete90/sun_w/si_liu/R/R-4.0.2/bin/R --vanilla CMD BATCH --quiet --no-save \
1a_saver_prepare_data.R 1a_saver_prepare_data.Rout


