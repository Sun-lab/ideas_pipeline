
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

if (length(args) < 5) {
  message("no enough arguments, using default values")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
  ncase    = 10      # case individuals
  nctrl    = 10      # control individuals
  ncell    = 120    # numbers of cells collected from each individuals.
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if(ncell == 0){
  UNEQ_N_CELL = TRUE
}else{
  UNEQ_N_CELL = FALSE
}

if(UNEQ_N_CELL){
  config = sprintf("ncase_%d_nctrl_%d_unequal_n_cell", ncase, nctrl)
}else{
  config = sprintf("ncase_%d_nctrl_%d_ncell_%d", ncase, nctrl, ncell)
}

config = sprintf("%s_fold_mean_%.1f_var_%.1f", config, r_mean, r_var)
config

# ---------------------------------------------------------------
# additional parameters
# ---------------------------------------------------------------

nCore = 6      # number of cores for multi-core computation
nall  = ncase + nctrl

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
library(transport)

library(data.table)
library(pryr)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

library(ideas)

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

# ---------------------------------------------------------------
# load data
# ---------------------------------------------------------------

sim_data     = readRDS(sprintf("data/sim_data_%s.rds", config))
count_matrix = sim_data$count_matrix
meta_cell    = sim_data$meta_cell
meta_ind     = sim_data$meta_ind
gene_index   = sim_data$gene_index

ls()
EE_index   = gene_index$EE_index
mean_index = gene_index$mean_index
var_index  = gene_index$var_index

dim(count_matrix)
count_matrix[1:3,1:6]

dim(meta_cell)
meta_cell[1:2,]

dim(meta_ind)
meta_ind[1:2,]

rm(sim_data)
gc()

# ---------------------------------------------------------------
# 1. DESeq2 analysis 
# ---------------------------------------------------------------
# We calculate count_matrix_bulk by adding up raw counts of
# all cells of an individual within a gene
# ---------------------------------------------------------------

# individual level info
u_ind = unique(meta_cell$individual)
length(u_ind)
u_ind[1:2]

count_matrix_bulk = matrix(nrow = nrow(count_matrix),
                         ncol = length(u_ind))
rownames(count_matrix_bulk) = rownames(count_matrix)
colnames(count_matrix_bulk) = u_ind

for (i_ind in 1:length(u_ind)) {
  cur_ind   = u_ind[i_ind]
  cur_ind_m = count_matrix[, meta_cell$individual == cur_ind]
  count_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
}

dim(count_matrix_bulk)
count_matrix_bulk[1:3,1:5]

meta_ind$phenotype = as.factor(meta_ind$phenotype)
dim(meta_ind)
meta_ind[1:2,]

# DESeq2 for all the genes
dds = DESeqDataSetFromMatrix(countData = count_matrix_bulk,
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

# ---------------------------------------------------------------
# 2. IDEAS 
# ---------------------------------------------------------------

var2test      = "phenotype"
var2adjust    = "RIN"
var2test_type = "binary"
var_per_cell  = "cell_rd"

dist_list = list()

for(fit_method in c("zinb", "kde")){
  for(d_metric in c("Was", "JSD")){
    message(sprintf("fit_method: %s, d_metric: %s\n", fit_method, d_metric))
    message(date())
    
    label = paste(fit_method, d_metric, sep="_")
    fnm = sprintf("data/dist_%s_%s.rds", label, config)
    if(file.exists(fnm)){
      dist_list[[label]] = readRDS(fnm)
    }else{
      dist1 = ideas_dist(count_matrix, meta_cell, meta_ind, 
                         var_per_cell, var2test, var2adjust, 
                         var2test_type, d_metric = d_metric, 
                         fit_method = fit_method)
      dist_list[[label]] = dist1
      saveRDS(dist1, fnm)
    }
  }
}

date()

lapply(dist_list, dim)

dist_list$zinb_Was[1,1:3,1:3]
dist_list$zinb_JSD[1,1:3,1:3]

# ---------------------------------------------------------------
# STEP 2: pval calculation 
# ---------------------------------------------------------------

n_gene = nrow(count_matrix)
y = as.numeric(meta_ind$phenotype==1)

pval_KR = matrix(NA, nrow=n_gene, ncol=length(dist_list))
rownames(pval_KR) = rownames(count_matrix)
colnames(pval_KR) = paste("KR", names(dist_list), sep="_")

for(k in 1:length(dist_list)){
  message(names(dist_list)[k])
  message(date())
  dist_k  = dist_list[[k]]
  pval_KR[,k] = foreach(i_g = 1:dim(dist_k)[1], .combine = "c") %dopar% {
    Di = dist_k[i_g,,]
    if(any(is.na(Di))){
      pval = NA 
    }else{
      Ki = D2K(Di)
      m1 = MiRKAT(y = y, X = meta_ind$RIN, Ks = Ki, out_type = "D", 
                  method = "permutation")
      pval = m1$p_values
    }
    pval
  }
}
date()

dim(pval_KR)
pval_KR[1:2,]


n_perm = 999
r.seed = 903

pval_PS = matrix(NA, nrow=n_gene, ncol=length(dist_list))
rownames(pval_PS) = rownames(count_matrix)
colnames(pval_PS) = names(dist_list)
colnames(pval_PS) = paste("PS", names(dist_list), sep="_")

for(k in 1:length(dist_list)){
  message(names(dist_list)[k])
  message(date())
  dist_k  = dist_list[[k]]
  pval_PS[,k] = permanova(dist_k, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm, 
                          r.seed=r.seed, residulize.x = FALSE)
}
date()

summary(pval_KR)
summary(pval_PS)

round(cor(-log10(pval_KR), use="pair"),2)
round(cor(-log10(pval_PS), use="pair"),2)

round(cor(-log10(pval_KR), -log10(pval_PS), use="pair"),2)

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

rds = colSums(count_matrix)
med_rds = median(rds)
summary(rds)
med_rds

dim(count_matrix)
count_matrix[1:3,1:6]
count_matrix = t(t(count_matrix)/rds)*med_rds
dim(count_matrix)
count_matrix[1:3,1:6]
summary(colSums(count_matrix))

count_matrix_log = log2(1 + count_matrix) #log transformed data

dim(count_matrix_log)
count_matrix_log[1:3, 1:4]
cell_id = colnames(count_matrix_log)   # get the cell id from the data
gene_id = rownames(count_matrix_log)   # get the gene id from the data

diagnosis = as.character(meta_cell$phenotype) #
diagnosis[diagnosis == 1] = "Case"
diagnosis[diagnosis == 0] = "Control"

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey  = cell_id)

sca = FromMatrix(count_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta_cell$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta_cell$individual)
colData(sca)$RIN = meta_ind$RIN[match(meta_cell$individual, 
                                      meta_ind$individual)]

colData(sca)

rm(count_matrix_log)
gc()

getOption("mc.cores")

date()
b0 = zlm(formula = ~ diagnosis + cngeneson + RIN, sca = sca, 
         parallel = TRUE)
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
# 4. ranksum analysis 
# ---------------------------------------------------------------
#
# the input of ranksum analysis can be matrix.
#
# input:
# (1) the expression count matrix,
#     with each column represents a cell and each row represents a gene.
# (2) the meta data, including cell individual information.

# ---------------------------------------------------------------
 
date()
ranksum_pval=apply(count_matrix,1,
                   function(x) wilcox.test(x[diagnosis=="Case"], 
                                           x[diagnosis=="Control"])$p.value)
date()

rm(count_matrix)

length(ranksum_pval)
ranksum_pval[1:4]

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

geneType = rep("EE", n_gene)
geneType[mean_index] = "meanDE"
geneType[var_index]  = "varDE"

df1 = data.frame(geneType, pval_KR, pval_PS, deseq2_pval, 
                 mast_pval_glm, mast_pval_glmer,ranksum_pval)
dim(df1)
df1[1:2,]


pdf(sprintf("figures/pvalue_hist_%s.pdf", config), 
    width = 9, height = 9)
par(mfrow = c(3, 3), mar=c(5,4,2,1), pty = "s", bty="n")
for(k in 2:ncol(df1)){
  plot.hist(df1[,k], idx_grp, names(df1)[k])
}
dev.off()

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------


fun1 <- function(x, alpha){table(x<=alpha)/length(x)}
apply(df1[,-1], 2, function(v){tapply(v, df1$geneType, fun1, alpha=0.05)})
apply(df1[,-1], 2, function(v){tapply(v, df1$geneType, fun1, alpha=0.01)})

write.table(df1, file=sprintf("results/pval_%s.txt", config), append=FALSE, 
            quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

sessionInfo()

mem_used()
gc()

q(save = "no")
