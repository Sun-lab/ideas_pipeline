#!/bin/bash

for threshold in 0.05
do
        R CMD BATCH --no-save --no-restore\
        "--args threshold='${threshold}'" 6e_summarize_results_qvalue.R 6e_summarize_results_qvalue_${threshold}.Rout
done
