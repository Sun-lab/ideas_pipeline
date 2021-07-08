#!/bin/bash
# sbatch -c 8 --mem 32G step1b_MAST_l4_per.sh
for type in $( cat cell_types_mast_l4.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_MAST_perm.R step1b_Rout/step1b_MAST_perm_${type}.Rout
done



