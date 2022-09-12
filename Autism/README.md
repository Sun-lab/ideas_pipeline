# ideas and other methods on Autism data

Run analysis using Autism data from paper Velmeshev et al.2019[[1]](#1). 

The goal is to do gene differential analysis in terms of case v.s. control.


DCA related preparation steps (this is needed for DCA\_direct):

[step1a_dca_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_dca_prepare_data.R)

[step1a_dca.sh](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_dca.sh)

**Note:** The output files from DCA for Autism data were got before 2021. DCA code had updates and as of Sept. 2021, there are two differences: 

(1) the command line to run DCA code for our setting changes to https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca.sh (for Autism data, PFC_all.csv should be used instead of all.csv)

(2) DCA output now longer directly provide mean_norm.tsv. We need to use the original count matrix and the mean.tsv output file from DCA to recover mean_norm.tsv. This step can be done using the code in https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca_recover_mean_norm.R.

[step1a_split_dca_outputs.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_split_dca_outputs.R) (note that in line 68, the pi.tsv file corresponds to dropout.tsv in the output files of DCA as of Sept.2021.)

SAVER related preparation steps (this is need if want to try saver\_direct):

[step1a_saver_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_saver_prepare_data.R)

### DESeq2

[step1b_DESeq2.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1b_DESeq2.R) DESeq2 without covariates and with covariates.

### Rank sum

[step1b_ranksum.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1b_ranksum.R).

### MAST

[step1b_MAST.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1b_MAST.R) MAST glm and glmer.


### Ideas\_nb

[step1c_ideas.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1c_ideas.R) runs MiRKAT and permanovas on distance matrix calculated based on nb(negative binominal) approach.

### DCA\_direct

[step1d_dca_direct.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1d_dca_direct.R)runs MiRKAT and permanovas on distance matrix calculated based on dca_direct(based on outputs from DCA[[2]](#2)) approach.

Note that here when using `fit_method = "dca_direct"` in function `ideas_dist`, the variable `"rd"` in `var_per_cell` is actually not involved in the distance matrix computing, unlike when using `fit_method = "nb"`. This input item is included here only to keep the code format consistent for different `fit_method` options. 
 
### SAVER\_direct

[step1d_saver.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1d_saver.R) uses SAVER[[3]](#3) as denosing method instead of DCA, for computing the distance matrix. 

### combine results, do GSEA and get q values

[step1e_combine_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1e_combine_results.R) combines p value results from step1b (results only from approach with covariates), step1c and step1e.

[step1f_gsea.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1f_gsea.R) Gene set enrichment analysis of results from step1e. 

[step1l_qvalues.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1l_qvalues.R) Computes corresponding q values. 


### summary steps:

[step2a_summarize_results.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step2a_summarize_results.R) summarizes GSEA results and comparison between results from different methods.

[step2b_summarize_results_perm.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step2b_summarize_results_perm.R) type-I error through permutation. 

[step2d_summarize_results_qvalue.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step2d_summarize_results_qvalue.R) overlap between significant genes from different methods under q-value cutoff. 

### DCA mean and pseudo dispersion exploration:

[step3a_DCA_formula_helper.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step3a_DCA_formula_helper.R) helper file getting indivdual level mean and variance for each gene, based on cell-level parameter estimates given by DCA.

[step3b_DCA_formula_covariates_all.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step3b_DCA_formula_covariates_all.R) gets p-values for association between log(mean) and mild/severe status, given covariates, by linear regression.

[step10b_DCA_pseudo_over_dispersion_covariates_all.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step10b_DCA_pseudo_over_dispersion_covariates_all.R) computes pseudo dispersion parameters and does regression.

[step10k_DCA_formula_four_groups_covariates_ranksum.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step10k_DCA_formula_four_groups_covariates_ranksum.R) ranksum tests on log(mean) (and log(pseudo theta)) regression p-values among groups separated by DESeq2 and DCA_direct q-value cutoffs.



## References
<a id="1">[1]</a> 
Velmeshev, Dmitry, et al. "Single-cell genomics identifies cell type–specific molecular changes in autism." Science 364.6441 (2019): 685-689.

<a id="2">[2]</a> 
Eraslan, Gökcen, et al. "Single-cell RNA-seq denoising using a deep count autoencoder." Nature communications 10.1 (2019): 1-14.

<a id="3">[3]</a> 
Huang, Mo, et al. "SAVER: gene expression recovery for single-cell RNA sequencing." Nature methods 15.7 (2018): 539-542.
