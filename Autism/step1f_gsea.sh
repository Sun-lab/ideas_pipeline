#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1f_gsea.R step1f_Rout/step1f_gsea_${type}.Rout
done

