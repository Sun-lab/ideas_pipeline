#!/bin/bash

for type in CD8+Tcells_1
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" 1b_DESeq2_mild_severe.R 1b_Rout/1b_DESeq2_mild_severe_${type}.Rout
done
