
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

sca2 = SummarizedExperiment(count_matrix, colData = meta_cell)
#sca2=sca2[1:50,(50*1:160)]
gc()
date()
zinb = zinbFit(sca2, K=2)
date()

sca_zinb = zinbwave(sca2, fitted_model = zinb, K = 2, epsilon=1000,
                    observationalWeights = TRUE)
date()
gc()

dds = DESeqDataSet(sca_zinb, design = ~ phenotype)
dds = DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
wv_pval = results(dds)$pvalue
names(wv_pval)=rownames(dds)

summary(wv_pval)
table(wv_pval[EE_index] < 0.05)
table(wv_pval[EE_index] < 0.05)/length(EE_index)


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


df1=read.table(file=sprintf("results/pval_ranksum_%s.txt", config), header=TRUE)

df1 = cbind(df1, data.frame(wv_pval))
#data.frame(geneType, pval_KR, pval_PS, deseq2_pval, 
#                 mast_pval_glm, mast_pval_glmer,wv_pval)
dim(df1)
df1[1:2,]


pdf(sprintf("figures/pvalue_hist_wv_%s.pdf", config), 
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

write.table(df1, file=sprintf("results/pval_wv_%s.txt", config), append=FALSE, 
            quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

sessionInfo()

#mem_used()
gc()

q(save = "no")
