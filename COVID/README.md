Run analysis using COVID data from paper Schulte-Schrepping et al. 2020 [[1]](#1). Data available at [https://beta.fastgenomics.org/home](https://beta.fastgenomics.org/home) by searching key word *Schulte-Schrepping* in data section. The data file used here is PBMC 10x data from cohort 1.  

The goal is to do gene differential analysis in terms of mild v.s. severe COVID diease status. 

Recover count data, fliter out genes appearing in less than 2000 cells and split into different cell types:

[0_separate_celltypes.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/0_separate_celltypes.R) 

DCA[[2]](#2) related preparation steps:

[1a_dca_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca_prepare_data.R)

[1a_dca.sh](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca.sh)

[1a_dca_recover_mean_norm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca_recover_mean_norm.R)

[1a_split_dca_outputs.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_split_dca_outputs.R)

Note on running DCA: 

(1) should install keras version 2.4 and tensorflow >=2.0 and <2.5, following the [requirements of DCA](https://github.com/theislab/dca/blob/master/setup.py),  not the new versions.  (2) current version (as of 09/26/2021) of DCA sets nb as default option. Need to set type to zinb-conddisp if zinb result is wanted. (3) current version does not output mean\_norm.tsv. Need to recover the mean\_norm.tsv file using mean.tsv file and the raw count file.  (4) current version does not provide meaningful column names for the dropout and dispersion matrix, need to borrow the column names (V1 followed by cell names) from mean\_norm matrix when splitting the matrices to celltypes. 

SAVER[[3]](#3) related preparation step:

[1a_saver_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_saver_prepare_data.R)

**DESeq2:**

[1b_DESeq2_mild_severe.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_DESeq2_mild_severe.R) 



**Rank sum test:**

[1b_ranksum.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_ranksum.R) 

Create permutated donor level labels:

[1b_permutation_label.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_permutation_label.R) Only involves 7 mild and 7 severe donors in each permutation to keep balance.

[1b_ranksum_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_ranksum_perm.R)


**MAST:**

[3b_MAST.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/3b_MAST.R) 


**IDEAS related:**

[5c_ideas.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5c_ideas.R)

[5d_dca_direct.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5d_dca_direct.R)

[5d_saver.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5d_saver.R)


**Summary steps:**

[5e_combine_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5e_combine_results.R) combines p-value results from different methods.

[5f_gsea.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5f_gsea.R) gene set enrichment analysis. 

[5l_qvalues.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5l_qvalues.R) gets q-values and number of genes significant under certain q-value cutoffs 0.1, 0.2, 0.3, and 0.4. 

[6a_summarize_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6a_summarize_results.R) prints out pathways significant under ajusted pvalue cutoff 0.05 from GSEA and does several other comparisons.

[6b_summarize_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6b_summarize_results.R) writes out detailed information on the significant pathways. 

[6c_summarize_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6c_summarize_results.R) fisher exact test among set of genes with p-value < 0.05 from different methods. 

[6e_summarize_results_qvalue.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6e_summarize_results_qvalue.R) gets number of genes with q-value < 0.05 from different methods and also gets the overlap counts and proportions(with respect to row name) between those from any pair of methods. 

**DCA mean and pseudo dispersion exploration:**

[7a_DCA_formula_helper.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/7a_DCA_formula_helper.R) helper file getting indivdual level mean and variance for each gene, based on cell-level parameter estimates given by DCA. 

[7b_DCA_formula_covariates_all.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/7b_DCA_formula_covariates_all.R) gets p-values for association between log(mean) and mild/severe status, given covariates, by linear regression. 

[7c_DCA_pseudo_over_dispersion_covariates_all.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/7c_DCA_pseudo_over_dispersion_covariates_all.R) computes pseudo dispersion parameters and does regression. 

[7d_DCA_formula_four_groups_covariates_ranksum.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/7d_DCA_formula_four_groups_covariates_ranksum.R) ranksum tests on log(mean) 
(and log(pseudo theta)) regression p-values among groups separated by DESeq2 and DCA\_direct q-value cutoffs.

**Mean expression level related exploration:**

[5f_neg_log_pvalue_vs_expression_level.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5f_neg_log_pvalue_vs_expression_level.R) gets mean expression level for each gene. 

[5f_gsea_expression.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5f_gsea_expression.R) GSEA when ranking genes by mean expression level. 

[7f_gsea_overlap_test.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/7f_gsea_overlap_test.R) counts and fisher exact test related quantities for the overlap between pathways significant under adjusted pvalue 0.05 from DESeq2, PS\_nb\_Was, PS\_dca\_direct\_Was and mean expression. 

**Permutations:**

[1b_permutation_label.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_permutation_label.R)

[1b_ranksum_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_ranksum_perm.R)

[3b_MAST_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/3b_MAST_perm.R)

[1b_DESeq2_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_DESeq2_perm.R)

[5c_ideas_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5c_ideas_perm.R)

[5d_dca_direct_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5d_dca_direct_perm.R)

[5d_saver_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/5d_saver_perm.R)

[6f_summarize_results_tpye_I.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6f_summarize_results_tpye_I.R)

[6g_summarize_results_FDR.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/6g_summarize_results_FDR.R)


## References
<a id="1">[1]</a> 
Schulte-Schrepping, Jonas, et al. "Severe COVID-19 is marked by a dysregulated myeloid cell compartment." Cell 182.6 (2020): 1419-1440.

<a id="2">[2]</a> 
Eraslan, Gökcen, et al. "Single-cell RNA-seq denoising using a deep count autoencoder." Nature communications 10.1 (2019): 1-14.

<a id="3">[3]</a> 
Huang, Mo, et al. "SAVER: gene expression recovery for single-cell RNA sequencing." Nature methods 15.7 (2018): 539-542.
