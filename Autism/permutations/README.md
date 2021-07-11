
[step11o_DESeq2_ten_p.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/permutations/step11o_DESeq2_ten_p.R) 10 permuations for DESeq2. 

[step11p_ideas_ten_p.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/permutations/step11p_ideas_ten_p.R) 10 permuations for ideas nb.

[step11q_dca_direct_ten_p.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/permutations/step11q_dca_direct_ten_p.R) 10 permutations for dca_direct.

[step11r_combine_results_ten_p.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/permutations/step11r_combine_results_ten_p.R) combines outputs from step11o, 11p, and 11q together, 
and writes the results to folder res](https://github.com/Sun-lab/ideas_pipeline/tree/main/Autism/permutations/res).

[step11t_hist_pvalues_ten_p.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/permutations/step11t_hist_pvalues_ten_p.R) draws histograms based on original results from 
step11e (the step name for the original runs was changed when reorganizing the code files, so the corresponding pvalue result files of step11e are also put under 
[res](https://github.com/Sun-lab/ideas_pipeline/tree/main/Autism/permutations/res) folder) and the combined premutation results from step11r. 
The histograms are put under folder [figure](https://github.com/Sun-lab/ideas_pipeline/tree/main/Autism/permutations/figures).
