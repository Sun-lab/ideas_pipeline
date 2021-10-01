
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

library(scDD)
library(SingleCellExperiment)


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
# 6. scDD analysis 
# ---------------------------------------------------------------
# We use the zinbwave to reduce the dimension, then apply DESeq2 for analysis
# all cells of an individual within a gene
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

cell_id = colnames(count_matrix)   # get the cell id from the data
gene_id = rownames(count_matrix)   # get the gene id from the data

gc()
date()
sca3 = SingleCellExperiment(list(normcounts=count_matrix))

colData(sca3)$condition = meta_cell$phenotype + 1
table(sca3$condition)

#sca3=sca3[1:100,50*(1:160)]
colData(sca3)
date()
sca_dd = scDD(sca3)
date()

res = results(sca_dd)
dim(res)
res[1:2,]

gc()

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

table(res$DDcategory, geneType)

pdf(sprintf("figures/pvalue_hist_scDD_%s.pdf", config), 
    width = 9, height = 9)
par(mfrow = c(3, 3), mar=c(5,4,2,1), pty = "s", bty="n")
plot.hist(res$nonzero.pvalue,  idx_grp, "nonzero.pvalue")
plot.hist(res$zero.pvalue,     idx_grp, "zero.pvalue")
plot.hist(res$combined.pvalue, idx_grp, "combined.pvalue")
dev.off()

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------

fun1 <- function(x, alpha){table(x<=alpha)/length(x)}
res_pval = res[,c("nonzero.pvalue", "zero.pvalue", "combined.pvalue")]
apply(res_pval, 2, function(v){tapply(v, geneType, fun1, alpha=0.05)})
apply(res_pval, 2, function(v){tapply(v, geneType, fun1, alpha=0.01)})

write.table(res, file=sprintf("results/res_scDD_%s.txt", config), 
            append=FALSE, quote=FALSE, sep="\t", row.names = FALSE, 
            col.names = TRUE)

sessionInfo()

mem_used()
gc()

q(save = "no")
