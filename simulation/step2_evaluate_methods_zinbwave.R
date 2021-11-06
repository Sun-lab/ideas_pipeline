
#########################################################################
#                                                                       #
#                                                                       #
#             Part II: Defferential Expression Analysis                 #
#                                                                       #
#                                                                       #
#########################################################################
# once simulated genes, we will do the Differential Expression Analysis
# here we implement the zinbwave

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) < 5) {
  message("no enough arguments, using default values")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
  ncase    = 10      # case individuals
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

library(zinbwave)
library(DESeq2)

library(data.table)
library(pryr)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

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
# 5. zinb_pval analysis 
# ---------------------------------------------------------------
# We use the zinbwave to estimate weights the dimension, then apply 
# DESeq2 given the weights
# ---------------------------------------------------------------

sca1 = SingleCellExperiment(list(counts=count_matrix), 
                            colData = meta_cell)
gc()
date()
sca_zinb = zinbwave(sca1, K = 2, observationalWeights = TRUE)
date()
gc()

dds = DESeqDataSet(sca_zinb, design = ~ phenotype)
dds = DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
res = results(dds)
dim(res)
res[1:2,]

summary(res$pvalue)
table(res$pvalue[EE_index] < 0.05)
table(res$pvalue[EE_index] < 0.05)/length(EE_index)

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------

geneType = rep("EE", nrow(count_matrix))
geneType[mean_index] = "meanDE"
geneType[var_index]  = "varDE"

fun1 <- function(x, alpha){table(x<=alpha)/length(x)}
tapply(res$pvalue, geneType, fun1, alpha=0.05)
tapply(res$pvalue, geneType, fun1, alpha=0.01)

res$gene = rownames(res)
res

write.table(res, file=sprintf("results/res_ZINB-WaVE_%s.txt", config), 
            append=FALSE, quote=FALSE, sep="\t", row.names = FALSE, 
            col.names = TRUE)

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

table(res$pvalue < 0.05, geneType)

pdf(sprintf("figures/pvalue_hist_zinbwave_%s.pdf", config), 
    width = 9, height = 3)
par(mfrow = c(1,3), mar=c(5,4,2,1), pty = "s", bty="n")
plot.hist(res$pvalue,  idx_grp, "ZINB-WaVE")
dev.off()


sessionInfo()

#mem_used()
gc()

q(save = "no")
