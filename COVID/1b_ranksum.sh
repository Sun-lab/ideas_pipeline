#!/bin/bash

for type in CD8+Tcells_1 NKcells
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" 1b_ranksum.R 1b_Rout/1b_ranksum_${type}.Rout
done



