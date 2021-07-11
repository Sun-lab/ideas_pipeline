#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1d_saver_perm.R step1d_Rout/step1d_saver_perm_${type}.Rout
done
