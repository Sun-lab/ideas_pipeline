#!/bin/bash


for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1e_combine_results.R step1e_Rout/step1e_combine_results_${type}.Rout
done

