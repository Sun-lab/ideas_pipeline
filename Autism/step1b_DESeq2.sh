#!/bin/bash

for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1b_DESeq2.R step1b_Rout/step1b_DESeq2_${type}.Rout
done



