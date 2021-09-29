Run analysis using COVID data from paper Schulte-Schrepping et al. 2020 [[1]](#1). Data available at [https://beta.fastgenomics.org/home](https://beta.fastgenomics.org/home) by searching key word *Schulte-Schrepping* in data section. The data file used here is PBMC 10x data from cohort 1. 

Recover count data, fliter out genes appearing in less than 2000 cells and split into different cell types:

[0_separate_celltypes.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/0_separate_celltypes.R) 

DCA[[2]](#2) related steps:

[1a_dca_prepare_data.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca_prepare_data.R)

[1a_dca.sh](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1a_dca.sh)

Note on running DCA: (1) should install keras version 2.4 and tensorflow >=2.0 and <2.5, following the [requirements of DCA](https://github.com/theislab/dca/blob/master/setup.py),  not the new versions. (2) current version (as of 09/26/2021) of DCA sets nb as default option. Need to set type to zinb-conddisp if zinb result is wanted. (3) current version does not output mean\_norm.tsv. 

DESeq2  :

[1b_DESeq2.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_DESeq2.R) 

[1b_DESeq2_mild_severe.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_DESeq2_mild_severe.R) 
Rank sum test:

[1b_ranksum.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_ranksum.R) 

MAST:

[1b_MAST.R](https://github.com/Sun-lab/ideas_pipeline/blob/main/COVID/1b_MAST.R) 







## References
<a id="1">[1]</a> 
Schulte-Schrepping, Jonas, et al. "Severe COVID-19 is marked by a dysregulated myeloid cell compartment." Cell 182.6 (2020): 1419-1440.

<a id="2">[2]</a> 
Eraslan, Gökcen, et al. "Single-cell RNA-seq denoising using a deep count autoencoder." Nature communications 10.1 (2019): 1-14.
