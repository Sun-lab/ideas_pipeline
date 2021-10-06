#!/bin/bash

for threshold in 0.01 0.05 0.10
do
        R CMD BATCH --no-save --no-restore\
        "--args threshold='${threshold}'" 4e_summarize_results_qvalue.R 4e_Rout/4e_summarize_results_qvalue_${threshold}.Rout
done
