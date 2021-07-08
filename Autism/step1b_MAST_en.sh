#!/bin/bash

for type in $( cat cell_types_mast_en.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_MAST_perm.R step1b_Rout/step1b_MAST_perm_${type}.Rout
done



