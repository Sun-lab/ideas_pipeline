# this file combines the pvalue results from step1b, step1c, 
# and step1d together to a dataframe

# files to load:
#   -- the 7 pvalue files
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

# DESeq2, MAST and ranksum pvalues did not filter out genes with 
# large proportion of zeros
# need to repeat the filtering out genes step to get the genes we want
pvals_deseq2     = fread(sprintf("res/step1b_DESeq2_%s_adj_covariates.tsv", grp))
pvals_nb_rd      = fread(sprintf("res/step1c_nb_pvals_%s.tsv", grp))
pvals_dca_direct = fread(sprintf("res/step1d_dca_direct_pvals_%s.tsv", grp))

pvals_mast_glm   = fread(sprintf("res/step1b_MAST_%s_glm.tsv", grp))
pvals_mast_glmer = fread(sprintf("res/step1b_MAST_%s_glmer.tsv", grp))

pvals_rank_sum   = fread(sprintf("res/step1b_ranksum_%s.tsv", grp))
pvals_saver      = fread(sprintf("res/step1d_saver_direct_pvals_%s.tsv", grp))

dim(pvals_deseq2)
dim(pvals_nb_rd)
dim(pvals_dca_direct)
dim(pvals_mast_glm)
dim(pvals_mast_glmer)
dim(pvals_rank_sum)
dim(pvals_saver)

pvals_deseq2[1:2,]
pvals_nb_rd[1:2,]
pvals_dca_direct[1:2,]
pvals_mast_glm[1:2,]
pvals_mast_glmer[1:2,]
pvals_rank_sum[1:2,]
pvals_saver[1:2,]

# add row name and gene column to help verify results after filtering
row.names(pvals_deseq2) = full_genes 
pvals_deseq2$gene = full_genes
pvals_deseq2 = pvals_deseq2[w2kp,]

pvals_mast_glm   = pvals_mast_glm[w2kp]
pvals_mast_glmer = pvals_mast_glmer[w2kp,]
pvals_rank_sum   = pvals_rank_sum[w2kp,]

stopifnot(all(pvals_deseq2$gene == pvals_nb_rd$gene))
stopifnot(all(pvals_deseq2$gene == pvals_dca_direct$gene))
stopifnot(all(pvals_deseq2$gene == pvals_mast_glm$V1))
stopifnot(all(pvals_deseq2$gene == pvals_mast_glmer$V1))
stopifnot(all(pvals_deseq2$gene == pvals_rank_sum$V1))
stopifnot(all(pvals_deseq2$gene == pvals_saver$gene))

pvals = cbind(pvals_nb_rd$gene, pvals_deseq2$pvalue, pvals_nb_rd[, 2:5], 
              pvals_dca_direct[, 2:5], pvals_mast_glm$V2, pvals_mast_glmer$V2,
              pvals_rank_sum$V2, pvals_saver[, 2:5])
names(pvals)[1] = "gene"
names(pvals)[2] = "DESeq2"
names(pvals)

names(pvals)[11] = "MAST_glm"
names(pvals)[12] = "MAST_glmer"
names(pvals)[13] = "rank_sum"

dim(pvals)
head(pvals)
# plot(-log10(pvals$PS_dca_direct_Was), -log10(pvals$PS_saver_direct_Was))

file.name = sprintf("res/step1e_pvals_%s.tsv", grp)
fwrite(pvals, file = file.name, sep = "\t")


gc()

sessionInfo()
q(save="no")
