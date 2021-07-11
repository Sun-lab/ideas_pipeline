#!/bin/bash
source /app/lmod/lmod/init/profile

#ml R
export OMP_NUM_THREADS=1

for type in $( cat cell_types.txt )
do
     sbatch -t 01-10:00:00 -c 12 \
        --mem=48000 -o rep11q.log \
        -p campus-new --job-name=step11q_$type \
        /fh/scratch/delete90/sun_w/si_liu/R/R-4.0.2/bin/R --vanilla CMD BATCH --quiet --no-save \
        "--args grp='${type}'" step11q_dca_direct_ten_p.R step11q_dca_direct_ten_p_${type}.Rout
done

