# step11r_combine_results_ten_p.R
# the main part of this file is carried over from 
# step11j_combine_results_p.R
#   - modify the output data frame to hold results from 
#     ten permutations



# the main part of code is carried over from step11e_combine_results.R
# except this time we summarize the results from permutation
# change the step name and add "_p" to files to save



# this file combines the pvalue results from step11g, step11h, 
# and step11i together to a dataframe

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

# all three results are in list format
pvals_deseq2     = readRDS(sprintf("res/step11o_DESeq2_%s_adj_covariates_ten_p.rds", 
                                 grp))
pvals_nb_rd      = readRDS(sprintf("res/step11p_nb_pvals_%s_ten_p.rds", grp))
pvals_dca_direct = readRDS(sprintf("res/step11q_dca_direct_pvals_%s_ten_p.rds", 
                                 grp))


length(pvals_deseq2)
length(pvals_nb_rd)
length(pvals_dca_direct)

# add row name and gene column to help verify results after filtering
test_genes = full_genes[w2kp]

DESeq2 = pvals_deseq2[[1]]$pvalue[w2kp]
pvals_nb_rd_cp = pvals_nb_rd[[1]][, 2:5]
pvals_dca_direct_cp = pvals_dca_direct[[1]][, 2:5]

for (j in 2:10){
  DESeq2 = c(DESeq2, pvals_deseq2[[j]]$pvalue[w2kp])
  pvals_nb_rd_cp = rbind(pvals_nb_rd_cp, pvals_nb_rd[[j]][, 2:5])
  pvals_dca_direct_cp = rbind(pvals_dca_direct_cp, pvals_dca_direct[[j]][, 2:5])
}

long_gene = rep(pvals_nb_rd[[1]]$gene, 10)

length(DESeq2)
dim(pvals_nb_rd_cp)
dim(pvals_dca_direct_cp)
length(long_gene)

mean(rep(test_genes, 10) == long_gene)

#pvals = cbind(pvals_nb_rd$gene, pvals_deseq2$pvalue, pvals_nb_rd[, 2:5], 
#              pvals_dca_direct[, 2:5])
pvals = cbind(long_gene, DESeq2, pvals_nb_rd_cp, 
              pvals_dca_direct_cp)

names(pvals)[1] = "gene"
names(pvals)[2] = "DESeq2"

dim(pvals)
head(pvals)
summary(pvals)

file.name = sprintf("res/step11r_pvals_%s_ten_p.tsv", 
                    grp)
fwrite(pvals, file = file.name, sep = "\t")