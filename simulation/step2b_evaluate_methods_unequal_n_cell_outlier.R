
#########################################################################
#                                                                       #
#                                                                       #
#             Part II: Defferential Expression Analysis                 #
#                                                                       #
#                                                                       #
#########################################################################
# once simulated genes, we will do the Differential Expression Analysis
# Here we will implement our method, plus the DESeq2 and the MAST analysis.


args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 2) {
  message("only take two arguments for r_mean and r_var")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

ncase    = 13      # case individuals
nctrl    = 10      # control individuals

config = sprintf("ncase_%d_nctrl_%d_unequal_n_cell_outlier", ncase, nctrl)
config = sprintf("%s_fold_mean_%.1f_var_%.1f", config, r_mean, r_var)
config

# ---------------------------------------------------------------
# additional parameters
# ---------------------------------------------------------------

nCore = 12           # number of cores for multi-core computation
nGeneMean  = 1000    # number of genes with different means in cases
nGeneVar   = 1000    # number of genes with different variance in cases
nGeneBlank = 6000   # number of genes equivalently expressed
nGeneTotal = nGeneMean + nGeneVar + nGeneBlank # total numbers of genes
nall       = ncase + nctrl


# ---------------------------------------------------------------
# initial setup
# ---------------------------------------------------------------

library(MASS)
library(emdbook)
library(moments)
library(MAST)
library(lme4)
library(DESeq2)
library(doParallel)
library(foreach)
library(doRNG)
library(MiRKAT)
library(reticulate)

library(data.table)
library(pryr)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

library(ideas)

# source("../functions/ZINB_fit_functions.R")
# source("../functions/kl_divergence_functions.R")
# source("../functions/Fstat_functions.R")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

# ---------------------------------------------------------------
# load data
# ---------------------------------------------------------------

gene_index = readRDS(sprintf("data/gene_index_%s.rds", config))
sim_matrix = readRDS(sprintf("data/sim_matrix_%s.rds", config))
meta_cell  = readRDS(sprintf("data/meta_%s.rds", config))
meta_ind   = readRDS(sprintf("data/meta_ind_%s.rds", config))

ls()
EE_index   = gene_index$EE_index
mean_index = gene_index$mean_index
var_index  = gene_index$var_index

dim(sim_matrix)
sim_matrix[1:2,1:5]

dim(meta_cell)
meta_cell[1:2,]

dim(meta_ind)
meta_ind[1:2,]

# ---------------------------------------------------------------
# 1. DESeq2 analysis 
# ---------------------------------------------------------------
# We calculate sim_matrix_bulk by adding up raw counts of
# all cells of an individual within a gene
# ---------------------------------------------------------------

# individual level info
u_ind = unique(meta_cell$individual)
length(u_ind)
u_ind[1:2]

sim_matrix_bulk = matrix(nrow = nrow(sim_matrix),
                         ncol = length(u_ind))
rownames(sim_matrix_bulk) = rownames(sim_matrix)
colnames(sim_matrix_bulk) = u_ind

for (i_ind in 1:length(u_ind)) {
  cur_ind   = u_ind[i_ind]
  cur_ind_m = sim_matrix[, meta_cell$individual == cur_ind]
  sim_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
}

meta_ind$phenotype = as.factor(meta_ind$phenotype)
dim(meta_ind)
meta_ind[1:2,]

# DESeq2 for all the genes
dds = DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                             colData = meta_ind,
                             design = ~ RIN + phenotype)

dds = DESeq(dds)
nms = resultsNames(dds)
nms

deseq2_pval = results(dds)$pvalue

rk = results(dds, name="RIN")
deseq2_pval_RIN = rk$pvalue

summary(deseq2_pval)
summary(deseq2_pval_RIN)
table(deseq2_pval[EE_index] < 0.05)
table(deseq2_pval[EE_index] < 0.05)/length(EE_index)

# DESeq2 only for the EE genes
dds0 = DESeqDataSetFromMatrix(countData = sim_matrix_bulk[EE_index,],
                              colData = meta_ind,
                              design = ~ RIN + phenotype)

dds0 = DESeq(dds0)
deseq2_pval0 = results(dds0)$pvalue

table(deseq2_pval0 < 0.05)
table(deseq2_pval0 < 0.05)/length(EE_index)

# ---------------------------------------------------------------
# 2. IDEAS 
# ---------------------------------------------------------------

count_matrix  = sim_matrix
var2test      = "phenotype"
var2adjust    = "RIN"
var2test_type = "binary"
var_per_cell  = "cell_rd"

date()
if(file.exists(sprintf("data/dist_zinb_%s.rds", config))){
  dist_zinb = readRDS(sprintf("data/dist_zinb_%s.rds", config))
}else{
  dist_zinb = ideas_dist(count_matrix, meta_cell, meta_ind, var_per_cell, 
                         var2test, var2adjust, var2test_type, 
                         fit_method = "zinb", per_cell_adjust = "NB")
  saveRDS(dist_zinb, sprintf("data/dist_zinb_%s.rds", config))
}
date()

dim(dist_zinb)
dist_zinb[1,1:3,1:3]

date()
if(file.exists(sprintf("data/dist_kde_%s.rds", config))){
  dist_kde = readRDS(sprintf("data/dist_kde_%s.rds", config))
}else{
  dist_kde = ideas_dist(count_matrix, meta_cell, meta_ind, var_per_cell, 
                        var2test, var2adjust, var2test_type, 
                        fit_method = "kde")
  saveRDS(dist_kde,  sprintf("data/dist_kde_%s.rds",  config))
}
date()

dim(dist_kde)
dist_kde[1,1:3,1:3]

# ---------------------------------------------------------------
# STEP 2: pval calculation 
# ---------------------------------------------------------------

y = as.numeric(meta_ind$phenotype==1)

date()
pval_M_zinb = foreach(i_g = 1:dim(dist_zinb)[1], .combine = "c") %dopar% {
  Di = dist_zinb[i_g,,]
  if(any(is.na(Di))){
    pval = NA 
  }else{
    Ki = D2K(Di)
    m1 = MiRKAT(y = y, X = meta_ind$RIN, Ks = Ki, out_type = "D", 
                method = "permutation")
    pval = m1$indivP
  }
  pval
}
date()

summary(pval_M_zinb)

date()
pval_M_kde = foreach(i_g = 1:dim(dist_kde)[1], .combine = "c") %dopar% {
  Di = dist_kde[i_g,,]
  if(any(is.na(Di))){
    pval = NA 
  }else{
    Ki = D2K(Di)
    m1 = MiRKAT(y = y, X = meta_ind$RIN, Ks = Ki, out_type = "D", 
                method = "permutation")
    pval = m1$indivP
  }
  pval
}
date()

summary(pval_M_kde)

n_perm = 999
r.seed = 819
delta  = 0.5

date()
pval_S1_zinb  = permanova(dist_zinb, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm, r.seed=r.seed, 
                          residulize.x = TRUE)
date()
pval_S0_zinb  = permanova(dist_zinb, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm, r.seed=r.seed, 
                          residulize.x = FALSE)
date()
pval_S1_kde   = permanova(dist_kde, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm, r.seed=r.seed, 
                          residulize.x = TRUE)
date()
pval_S0_kde   = permanova(dist_kde, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm, r.seed=r.seed, 
                          residulize.x = FALSE)
date()

table(pval_M_zinb < 0.01, pval_M_kde < 0.01)
table(pval_S0_zinb < 0.01, pval_S0_kde < 0.01)
table(pval_S1_zinb < 0.01, pval_S1_kde < 0.01)

table(pval_M_zinb < 0.01, pval_S0_zinb < 0.01)
table(pval_M_zinb < 0.01, pval_S1_zinb < 0.01)

table(pval_M_kde < 0.01, pval_S0_kde < 0.01)
table(pval_M_kde < 0.01, pval_S1_kde < 0.01)

# ---------------------------------------------------------------
# 3. MAST analysis 
# ---------------------------------------------------------------
#
# the input of MAST analysis can be matrix or SingleCellAssay.
#
# input:
# (1) the log-transformed expression count matrix,
#     with each column represents a cell and each row represents a gene.
# (2) the meta data, including cell and individual information.

# we get the p-values based on the Hurdle model ("H" model)
# ---------------------------------------------------------------

sim_matrix_log = log2(1 + sim_matrix) #log transformed data

dim(sim_matrix_log)
sim_matrix_log[1:3, 1:4]
cell_id = colnames(sim_matrix_log)   # get the cell id from the data
gene_id = rownames(sim_matrix_log)   # get the gene id from the data

diagnosis = as.character(meta_cell$phenotype) #
diagnosis[diagnosis == 1] = "Case"
diagnosis[diagnosis == 0] = "Control"

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey  = cell_id)

sca = FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta_cell$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta_cell$individual)
colData(sca)$RIN = meta_ind$RIN[match(meta_cell$individual, meta_ind$individual)]

colData(sca)

getOption("mc.cores")

date()
b0 = zlm(formula = ~ diagnosis + cngeneson + RIN, sca = sca, parallel = TRUE)
date()
b1 = zlm(formula = ~ diagnosis + (1 | ind) + cngeneson + RIN, sca = sca, 
         method = 'glmer', ebayes = FALSE, parallel = TRUE)
date()

b0
b1

date()
lrt0 = lrTest(b0, "diagnosis")
date()
lrt1 = lrTest(b1, "diagnosis")
date()

dim(lrt0)
lrt0[1,,]

dim(lrt1)
lrt1[1,,]

mast_pval_glm = apply(lrt0, 1, function(x){x[3,3]})
length(mast_pval_glm)
mast_pval_glm[1:4]

mast_pval_glmer = apply(lrt1, 1, function(x){x[3,3]})
length(mast_pval_glmer)
mast_pval_glmer[1:4]

# ---------------------------------------------------------------
# check p-value
# ---------------------------------------------------------------

idx_grp = list(meanDE=mean_index, varDE=var_index, EE=EE_index)

plot.hist <- function(pvals, idx_grp, label){
  for(k in 1:length(idx_grp)){
    idx  = idx_grp[[k]]
    main = paste(label, names(idx_grp)[k], sep=", ")
    hist(pvals[idx], main=main, xlab="p-value", breaks=50)
  }
}

pdf(sprintf("figures/pvalue_hist_%s.pdf", config), 
    width = 9, height = 9)
par(mfrow = c(3, 3), mar=c(5,4,2,1), pty = "s", bty="n")
plot.hist(pval_M_zinb,    idx_grp, "IDEAS-M-ZINB")
plot.hist(pval_M_kde,     idx_grp, "IDEAS-M-KDE")
plot.hist(pval_S0_zinb,   idx_grp, "IDEAS-S0-ZINB")
plot.hist(pval_S0_kde,    idx_grp, "IDEAS-S0-KDE")
plot.hist(pval_S1_zinb,   idx_grp, "IDEAS-S1-ZINB")
plot.hist(pval_S1_kde,    idx_grp, "IDEAS-S1-KDE")

plot.hist(deseq2_pval,     idx_grp, "DESeq2")
plot.hist(mast_pval_glm,   idx_grp, "MAST bayesglm")
plot.hist(mast_pval_glmer, idx_grp, "MAST glmer")
dev.off()

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------

geneType = rep("EE", nGeneTotal)
geneType[mean_index] = "meanDE"
geneType[var_index]  = "varDE"

df1 = data.frame(geneType, pval_M_zinb, pval_M_kde, pval_S0_zinb, 
                 pval_S0_kde, pval_S1_zinb, pval_S1_kde,
                 deseq2_pval, mast_pval_glm, mast_pval_glmer)
dim(df1)
df1[1:2,]

fun1 <- function(x, alpha){table(x<=alpha)/length(x)}
apply(df1[,-1], 2, function(v){tapply(v, df1$geneType, fun1, alpha=0.05)})

apply(df1[,-1], 2, function(v){tapply(v, df1$geneType, fun1, alpha=0.01)})

write.table(df1, file=sprintf("results/pval_%s.txt", config), append=FALSE, 
            quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

sessionInfo()

mem_used()
gc()

q(save = "no")
