#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1d_saver.R step1d_Rout/step1d_saver_${type}.Rout
done
