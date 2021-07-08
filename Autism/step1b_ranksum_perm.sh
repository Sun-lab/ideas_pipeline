#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_ranksum_perm.R step1b_Rout/step1b_ranksum_perm_${type}.Rout
done



