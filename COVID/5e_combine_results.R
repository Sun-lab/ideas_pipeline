# compared to 3e, this version combines results from 1b DESeq2, 1b ranksum test, 
# 3b MAST, 5c, 5d, allowing higher permutation number of times

# compared to 1e, this version combines results from 1b DESeq2, 1b ranksum test, 
# 3b MAST, 3c, 3d, considers more genes
# by setting filtering criterion to keep genes
# appearning in at last 10% of the cells



# this file combines the pvalue results from 1b, 1c, 
# and 1d together to a dataframe

# files to load:
#   -- the 7 pvalue files
#   -- the orignal full count and cell info matrices to filter 
#      out genes with too many zeros

# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'CD8+Tcells_1' as default.\n")
  grp = "CD8+Tcells_1"
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


data.dir = "../../ideas_data/COVID/PBMC_10x"



# DESeq2 and ranksum pvalues did not filter out genes with 
# large proportion of zeros (MAST for COVID has done that)
# need to repeat the filtering out genes step to get the genes we want
pvals_deseq2     = fread(sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_mild_severe.tsv", grp))
pvals_nb_rd      = fread(sprintf("res/5c_nb_pvals_%s.tsv", grp))
pvals_dca_direct = fread(sprintf("res/5d_dca_direct_pvals_%s.tsv", grp))

pvals_mast_glm   = fread(sprintf("res/3b_MAST_%s_glm.tsv", grp))
pvals_mast_glmer = fread(sprintf("res/3b_MAST_%s_glmer.tsv", grp))

pvals_rank_sum   = fread(sprintf("res/1b_ranksum_%s.tsv", grp))
pvals_saver      = fread(sprintf("res/5d_saver_direct_pvals_%s.tsv", grp))

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
dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
full_genes = rownames(dat)  
row.names(pvals_deseq2) = full_genes 
pvals_deseq2$gene = full_genes



# get filtering condition
cell_info = fread(file.path(data.dir, "meta.tsv"), 
                  stringsAsFactors=TRUE)

meta = cell_info[match(colnames(dat), cell_info$cell),]
meta_covid = meta[which(meta$group_per_sample != "control"),]

df_donor = as.data.frame(table(meta_covid$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 10)]

meta2kp = meta_covid[which(meta_covid$donor %in% donor2kp),]
cell2kp_index = which(meta$cell %in% meta2kp$cell)

dat1 = dat[, cell2kp_index]

n.zeros = rowSums(dat1 == 0)
summary(n.zeros)

w2kp = which(n.zeros < 0.9*ncol(dat1))
length(w2kp)


# subset pvalues for genes to keep
pvals_deseq2 = pvals_deseq2[w2kp,]
pvals_rank_sum   = pvals_rank_sum[w2kp,]

table(pvals_deseq2$gene == pvals_nb_rd$gene)
table(pvals_deseq2$gene == pvals_dca_direct$gene)
table(pvals_deseq2$gene == pvals_mast_glm$V1)
table(pvals_deseq2$gene == pvals_mast_glmer$V1)
table(pvals_deseq2$gene == pvals_rank_sum$V1)
table(pvals_deseq2$gene == pvals_saver$gene)

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

file.name = sprintf("res/5e_pvals_%s.tsv", grp)
fwrite(pvals, file = file.name, sep = "\t")


gc()

sessionInfo()
q(save="no")
