# ideas methods on Autism data

### step1a_dca_prepare_data.R

Prepare input data for DCA. Only needed if want to try dca_direct in step1d.

### step1a_dca.sh

Command line for running DCA. Only needed if want to try dca_direct in step1d.

### step1a_split_dca_outputs.R

Splitting the outputs of dca according to difference cell types, to prepare part of inputs for step1d. Only needed if want to try dca_direct in step1d. 

### step1b_DESeq2.R

DESeq2 without covariates and with covariates.

### step1c_ideas.R

Run MiRKAT and permanovas on distance matrix calculated based on nb(negative binominal) approach.

### step1d_dca_direct.R

Run MiRKAT and permanovas on distance matrix calculated based on dca_direct(based on outputs from DCA) approach.

Note that here when using `fit_method = "dca_direct"` in function `ideas_dist`, the variable `"rd"` in `var_per_cell` is actually not involved in the distance matrix computing, unlike when using `fit_method = "nb"`. This input item is included here only to keep the code format consistent for different `fit_method` options. 
 
### step1e_combine_results.R

Combine p value results from step1b (results only from approach with covariates), step1c and step1e.

### step1f_gsea.R

Gene set enrichment analysis of results from step1e. 

### step1l_qvalues.R

Compute corresponding q values. 

### cell_types.txt

List all the 17 cell types for looping through using .sh files when needed.
