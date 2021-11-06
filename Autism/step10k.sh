#!/bin/bash

R CMD BATCH "--args deseq2_qcut=0.2 dca_direct_qcut=0.1" \
            step10k_DCA_formula_four_groups_covariates_ranksum.R \
            step10k_DCA_formula_four_groups_covariates_ranksum_0.2_0.1.Rout

R CMD BATCH "--args deseq2_qcut=0.2 dca_direct_qcut=0.2" \
            step10k_DCA_formula_four_groups_covariates_ranksum.R \
            step10k_DCA_formula_four_groups_covariates_ranksum_0.2_0.2.Rout
