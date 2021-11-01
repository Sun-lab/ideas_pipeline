
# illustrate a few examples

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) < 5) {
  message("no enough arguments, using default values")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
  ncase    = 5       # case individuals
  nctrl    = 5       # control individuals
  ncell    = 360     # numbers of cells collected from each individuals.
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
library(stringr)

library(data.table)
library(pryr)
library(ggplot2)
library(ggpubr)
theme_set(theme_classic())

library(ideas)
library(reshape2)

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
# read in p-values
# ---------------------------------------------------------------

pvals = fread(sprintf("results/pval_%s.txt", config))
dim(pvals)
pvals[1:2,]

table(pvals$geneType)

pvals_rank_sum = fread(sprintf("results/pval_ranksum_%s.txt", config))
dim(pvals_rank_sum)
pvals_rank_sum[1:2,]

stopifnot(max(abs(pvals$PS_zinb_Was - 
                    pvals_rank_sum$PS_zinb_Was), na.rm=TRUE) < 1e-10)
pvals = pvals_rank_sum

# ---------------------------------------------------------------
# read in scDD p-values
# ---------------------------------------------------------------

pvals_scDD = fread(sprintf("results/res_scDD_%s.txt", config))
dim(pvals_scDD)
pvals_scDD[1:2,]

pvals$scDD = pvals_scDD$combined.pvalue
pvals$scDD_category = pvals_scDD$DDcategory

t1 = table(pvals$geneType, pvals$scDD_category)
t1

scDD_cat = melt(t1, varnames = c("gene_type", "scDD_category"))

g1 = ggplot(scDD_cat %>% group_by(gene_type) %>% 
         mutate(rel_freq = round(value/sum(value),2)), 
       aes(x = gene_type, y = rel_freq, 
           fill = scDD_category, cumulative = TRUE)) +
  geom_col() +
  geom_text(aes(label = paste0(rel_freq*100,"%")), 
            position = position_stack(vjust = 0.5))

pdf(sprintf("figures/scDD_cat_%s.pdf", config), width=4, height=4)
print(g1)
dev.off()

# ---------------------------------------------------------------
# read in ZINB-WaVE p-values
# ---------------------------------------------------------------

pvals_zw = fread(sprintf("results/res_ZINB-WaVE_%s.txt", config))
dim(pvals_zw)
pvals_zw[1:2,]

pvals$ZINB_WaVE = pvals_zw$pvalue

# ---------------------------------------------------------------
# find a few examples to illustrate
# ---------------------------------------------------------------

w2use = which(pvals$geneType == "varDE" & pvals$deseq2_pval > 0.1 & 
                pvals$PS_zinb_Was <= 0.001)
length(w2use)
w2use
pvals[w2use,]

if(length(w2use) >= 2){
  row_id = w2use[2]
  row_id
  
  # ------------------------------------------------------------------------
  # extract gene expression in cells and bulk samples
  # ------------------------------------------------------------------------
  
  ct_cell = count_matrix[row_id,]
  length(ct_cell)
  print(table(ct_cell))
  meta_ind$diagnosis = as.factor(meta_ind$phenotype)
  
  ct_ind = tapply(ct_cell, as.character(meta_cell$individual), sum)
  mat1 = match(meta_ind$individual, names(ct_ind))
  table(names(ct_ind)[mat1] == meta_ind$individual)
  meta_ind[["gene1"]] = as.numeric(ct_ind)[mat1]
  
  nb2 = glm.nb(gene1 ~ RIN + diagnosis, data=meta_ind)
  print(summary(nb2))
  
  # ------------------------------------------------------------------------
  # boxplot of bulk gene expression vs. diagnosis
  # ------------------------------------------------------------------------
  
  p1 = ggplot(meta_ind, aes(x=diagnosis, y=log10(gene1+0.5), col=diagnosis)) + 
    geom_boxplot() + labs(y="log10(ind. level counts)")
  p1 = p1 + geom_jitter(shape=16, position=position_jitter(0.2))
  
  # ------------------------------------------------------------------------
  # density plot of cell level gene expression vs. diagnosis
  # ------------------------------------------------------------------------
  
  mat2 = match(meta_cell$individual, meta_ind$individual)
  meta_cell$diagnosis = as.factor(meta_ind$diagnosis[mat2])
  
  df_test = meta_cell
  df_test$count = ct_cell
  mat2 = match(meta_cell$individual, meta_ind$individual)
  df_test$diagnosis = as.factor(df_test$diagnosis)
  
  table(df_test$count)
  tb0 = table(df_test$count[which(df_test$phenotype==0)])
  tb1 = table(df_test$count[which(df_test$phenotype==1)])
  
  tb0
  tb1
  
  df_test$count[which(df_test$count >= 7)] = 7
  
  p3 = ggplot(df_test, aes(x=count, col=diagnosis, line_type=individual)) + 
    geom_freqpoly(binwidth=1, closed="left") + guides(color = FALSE) + 
    xlim(0, 7) + ylab("frequency")
  
  gg0 = ggarrange(p1, p3, ncol=2, nrow=1, widths = c(1.5, 2))
  
  pdf(sprintf("figures/ex1_%s.pdf", config), width=6.5, height=2.5)
  print(gg0)
  dev.off()
}

# ------------------------------------------------------------------------
# extract gene expression in cells and bulk samples
# ------------------------------------------------------------------------

cal.power <- function(x, geneType){
  tapply(x, geneType, function(x){sum(x < 0.05, na.rm=TRUE)/sum(!is.na(x))})
}

ms = c("PS_zinb_Was", "PS_kde_Was", "deseq2_pval", "mast_pval_glm", 
       "mast_pval_glmer", "ranksum_pval", "scDD", "ZINB_WaVE")
powers = apply(pvals[,..ms], 2, cal.power, geneType=pvals$geneType)

print(config)
print(powers)

gg = melt(powers)

names(gg) = c("geneType", "method", "power")
gg$method = gsub("deseq2_pval", "DEseq2", gg$method)
gg$method = gsub("PS_zinb_Was", "IDEAS_ZINB",  gg$method)
gg$method = gsub("PS_kde_Was",  "IDEAS_KDE", gg$method)
gg$method = gsub("mast_pval_glmer", "MAST_glmer", gg$method)
gg$method = gsub("mast_pval_glm", "MAST", gg$method)
gg$method = gsub("ranksum_pval", "Rank-sum", gg$method)
table(gg$method)
gg$method = factor(gg$method, 
                   levels = c("Rank-sum", "MAST", "scDD", "ZINB_WaVE",  
                              "MAST_glmer", "DEseq2", "IDEAS_ZINB", "IDEAS_KDE"))

g1 = ggplot(subset(gg, geneType %in% c("EE")), 
            aes(x=geneType, y=power, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Paired") + 
  geom_hline(yintercept=0.05, col="red") + 
  ylab("type I error")

g2 = ggplot(subset(gg, geneType %in% c("meanDE", "varDE")), 
            aes(x=geneType, y=power, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Paired") + 
  geom_hline(yintercept=0.05, col="red")

gg1 = ggarrange(g1, g2, ncol = 2, nrow = 1, widths=c(1.25,2), 
                common.legend = TRUE, legend = "top")

pdf(sprintf("figures/power_%s.pdf", config), width=5.2, height=2.5)
print(gg1)
dev.off()

sessionInfo()

mem_used()
gc()

q(save = "no")
