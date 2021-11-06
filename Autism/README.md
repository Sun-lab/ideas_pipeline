# ideas and other methods on Autism data

DCA related preparation steps (this is needed for DCA\_direct):

[step1a_dca_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_dca_prepare_data.R)

[step1a_dca.sh](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_dca.sh)

[step1a_split_dca_outputs.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1a_split_dca_outputs.R)

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

[step1d_dca_direct.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1d_dca_direct.R)runs MiRKAT and permanovas on distance matrix calculated based on dca_direct(based on outputs from DCA) approach.

Note that here when using `fit_method = "dca_direct"` in function `ideas_dist`, the variable `"rd"` in `var_per_cell` is actually not involved in the distance matrix computing, unlike when using `fit_method = "nb"`. This input item is included here only to keep the code format consistent for different `fit_method` options. 
 
### SAVER\_direct

[step1d_saver.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1d_saver.R) uses SAVER as denosing method instead of DCA, for computing the distance matrix. 

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





