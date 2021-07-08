# this file combines the pvalue results from step1b, step1c, 
# and step1d together to a dataframe

# files to load:
#   -- the 3 pvalue files
#   -- the orignal full count matrix to filter out genes with 
#      too many zeros

# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
}else{
  eval(parse(text=args[[1]]))
}

grp

# ========================================================================
# libraries and path
# ========================================================================

library(MASS)
library(Matrix)
library(data.table)


data.dir  = "./data"

full_dat1 = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
full_genes = rownames(full_dat1)  

n.zeros = rowSums(full_dat1 == 0)
summary(n.zeros)
w2kp = which(n.zeros < 0.8*ncol(full_dat1))
length(w2kp)

# DESeq2 pvalues did not filter out genes with large proportion of zeros
# need to repeat the filtering out genes step to get the genes we want
pvals_deseq2     = fread(sprintf("res/step1b_DESeq2_%s_adj_covariates.tsv", grp))
pvals_nb_rd      = fread(sprintf("res/step1c_nb_pvals_%s.tsv", grp))
pvals_dca_direct = fread(sprintf("res/step1d_dca_direct_pvals_%s.tsv", grp))


# add row name and gene column to help verify results after filtering
row.names(pvals_deseq2) = full_genes 
pvals_deseq2$gene = full_genes
pvals_deseq2 = pvals_deseq2[w2kp]


pvals = cbind(pvals_nb_rd$gene, pvals_deseq2$pvalue, pvals_nb_rd[, 2:5], 
              pvals_dca_direct[, 2:5])
names(pvals)[1] = "gene"
names(pvals)[2] = "DESeq2"

dim(pvals)
head(pvals)

file.name = sprintf("res/step1e_pvals_%s.tsv", 
                    grp)
fwrite(pvals, file = file.name, sep = "\t")