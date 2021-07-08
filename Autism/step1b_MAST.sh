#!/bin/bash

for type in $( cat cell_types_mast.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_MAST.R step1b_Rout/step1b_MAST_${type}.Rout
done



