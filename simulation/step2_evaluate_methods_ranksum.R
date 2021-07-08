
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
  ncase    = 13      # case individuals
  nctrl    = 10      # control individuals
  ncell    = 360    # numbers of cells collected from each individuals.
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

library(doParallel)
library(foreach)
library(doRNG)

library(data.table)
library(ggplot2)
theme_set(theme_bw())

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

rm(sim_data)
gc()


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

diagnosis = as.character(meta_cell$phenotype) #
diagnosis[diagnosis == 1] = "Case"
diagnosis[diagnosis == 0] = "Control"
 
date()
ranksum_pval=apply(count_matrix,1,function(x) wilcox.test(x[diagnosis=="Case"],x[diagnosis=="Control"])$p.value)
date()

#rm(count_matrix)

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

geneType = rep("EE", nrow(count_matrix))
geneType[mean_index] = "meanDE"
geneType[var_index]  = "varDE"


df1=read.table(file=sprintf("results/pval_%s.txt", config), header=TRUE)

df1 = cbind(df1, data.frame(ranksum_pval))
#data.frame(geneType, pval_KR, pval_PS, deseq2_pval, 
#                 mast_pval_glm, mast_pval_glmer,ranksum_pval)
dim(df1)
df1[1:2,]


pdf(sprintf("figures/pvalue_hist_ranksum_%s.pdf", config), 
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

write.table(df1, file=sprintf("results/pval_ranksum_%s.txt", config), append=FALSE, 
            quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

sessionInfo()

#mem_used()
gc()

q(save = "no")
