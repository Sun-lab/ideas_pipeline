#!/bin/bash

R CMD BATCH 3d_dca_direct.R

R CMD BATCH 3d_saver.R

R CMD BATCH 3e_combine_results.R

R CMD BATCH 3f_gsea.R

R CMD BATCH 3l_qvalues.R

R CMD BATCH 4a_summarize_results.R


for threshold in 0.01 0.05 0.10
do
        R CMD BATCH --no-save --no-restore\
        "--args threshold='${threshold}'" 4e_summarize_results_qvalue.R 4e_Rout/4e_summarize_results_qvalue_${threshold}.Rout
done
