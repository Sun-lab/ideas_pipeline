#!/bin/bash

R CMD BATCH "--args deseq2_qcut=0.05 dca_direct_qcut=0.05" \
            7d_DCA_formula_four_groups_covariates_ranksum.R \
            7d_log_formula_four_groups_0.05_0.05_covariates_ranksum.Rout

R CMD BATCH "--args deseq2_qcut=0.1 dca_direct_qcut=0.05" \
            7d_DCA_formula_four_groups_covariates_ranksum.R \
            7d_log_formula_four_groups_0.10_0.05_covariates_ranksum.Rout

R CMD BATCH "--args deseq2_qcut=0.2 dca_direct_qcut=0.1" \
            7d_DCA_formula_four_groups_covariates_ranksum.R \
            7d_log_formula_four_groups_0.20_0.10_covariates_ranksum.Rout

R CMD BATCH "--args deseq2_qcut=0.2 dca_direct_qcut=0.2" \
            7d_DCA_formula_four_groups_covariates_ranksum.R \
            7d_log_formula_four_groups_0.20_0.20_covariates_ranksum.Rout
