#!/bin/bash
for type in $( cat cell_types.txt )
do
        R CMD BATCH --no-save --no-restore\
        "--args grp='${type}'" step1c_ideas.R step1c_Rout/step1c_ideas_${type}.Rout
done
