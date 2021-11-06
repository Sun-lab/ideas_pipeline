#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1d_dca_direct.R step1d_Rout/step1d_dca_direct_${type}.Rout
done
