# This code is modified from
# step3e_DCA_formula_pmf_three_groups_covariates_ranksum.R
# to separate DESeq2 q value and dca_direct q value 
# at two pairs of cutoffs:
# (1) DESeq2 q value 0.2, dca_direct q value 0.2
# (2) DESeq2 q value 0.2, dca_direct q value 0.1
# these two are taken as inputs to this file


## current version handles L2_3 only


args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 2) {
  message("two arguments are expected, use 0.2 and 0.1 as default.\n")
  deseq2_qcut = 0.2
  dca_direct_qcut = 0.1
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

deseq2_qcut = as.numeric(deseq2_qcut)
deseq2_qcut

dca_direct_qcut = as.numeric(dca_direct_qcut)
dca_direct_qcut



library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)

library(grid)
library(gridExtra)

theme_set(theme_classic())





data.dir  = "./data"

grp = "PFC_L2_3"
grp1 = "L2_3"



## data processing

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"), 
                  stringsAsFactors=TRUE)
dim(cell_info)
cell_info[1:2,]

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------


dat1 = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
full_genes = rownames(dat1)  



dim(dat1)
dat1[1:5,1:4]


# ------------------------------------------------------------------------
# filter out genes with too many zero's
# ------------------------------------------------------------------------

n.zeros = rowSums(dat1 == 0)
summary(n.zeros)

0.6*ncol(dat1)
0.8*ncol(dat1)

table(n.zeros < 0.6*ncol(dat1))
table(n.zeros < 0.8*ncol(dat1))
table(n.zeros < 0.9*ncol(dat1))

w2kp = which(n.zeros < 0.8*ncol(dat1))
dat1 = dat1[w2kp,]

dim(dat1)
dat1[1:5,1:4]



# -------------------------------------------------------------------
# load means and variances pvalues from log transformation
# and linear regression on covariates
# from formula and pmf based approach
# -------------------------------------------------------------------

dca_pvalues = read.csv(
  "res/step3b_formula_covariates_pvals_all.csv", 
  header = TRUE)
dim(dca_pvalues)
dca_pvalues[1:2, ]

theta_pvalues = read.csv(
  "res/step10b_formula_covariates_pvals_all.csv", 
  header = TRUE)
dim(theta_pvalues)
theta_pvalues[1:2, ]

dca_pvalues$theta_formula = theta_pvalues$theta_formula
dim(dca_pvalues)
dca_pvalues[1:2, ]

# -------------------------------------------------------------------
# read in q value information
# prepare for comparing gene groups
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

i = 8
ct1 = cell_types[i]

methods = c("DESeq2", "rank_sum", "MAST_glm", "MAST_glmer", 
            "PS_nb_Was", "PS_dca_direct_Was", "PS_saver_direct_Was")

q1  =  fread(sprintf("res/step1l_qvals_%s.tsv", ct1))
dim(q1)
table(q1$gene == rownames(dat1))

threshold_1 = 0.1
threshold_2 = 0.2
threshold_3 = 0.3 

check_length_1 <- function(v){
  length(which(v < threshold_1))
}

check_length_2 <- function(v){
  length(which(v < threshold_2))
}

check_length_3 <- function(v){
  length(which(v < threshold_3))
}

apply(q1[, -1], 2, check_length_1)

apply(q1[, -1], 2, check_length_2)

apply(q1[, -1], 2, check_length_3)



# ------------------------------------------------------------
# pvalues on log(mean) and log(variance) from dca formula 
# ------------------------------------------------------------

  
deseq2_index = q1$DESeq2 < deseq2_qcut

dca_direct_index = q1$PS_dca_direct_Was < dca_direct_qcut


group1_index = deseq2_index & dca_direct_index
group2_index = (!deseq2_index) & dca_direct_index
group3_index = deseq2_index & (!dca_direct_index)
group4_index = (!deseq2_index) & (!dca_direct_index)

sum(group1_index)
sum(group2_index)
sum(group3_index)
sum(group4_index)

sum(sum(group1_index) + sum(group2_index) +
    sum(group3_index) + sum(group4_index))

group1_genes = q1$gene[group1_index]
group1_genes[1:10]

group2_genes = q1$gene[group2_index]
group2_genes[1:10]

group3_genes = q1$gene[group3_index]
group3_genes[1:10]

group4_genes = q1$gene[group4_index]
group4_genes[1:10]

rows1 = match(group1_genes, row.names(dat1))
rows2 = match(group2_genes, row.names(dat1))
rows3 = match(group3_genes, row.names(dat1))  
rows4 = match(group4_genes, row.names(dat1))  

table(which(group1_index==TRUE) == rows1)
table(which(group2_index==TRUE) == rows2)
table(which(group3_index==TRUE) == rows3)
table(which(group4_index==TRUE) == rows4)

  
# for each gene in interest
# get log(mean) and log(theta) pvalue for 
# each individual in each group
# test the difference in mean between case and control
#                     '''theta '''

group1_mean_pval = dca_pvalues$mean_formula[rows1]
group1_theta_pval = dca_pvalues$theta_formula[rows1]

group2_mean_pval = dca_pvalues$mean_formula[rows2]
group2_theta_pval = dca_pvalues$theta_formula[rows2] 

group3_mean_pval = dca_pvalues$mean_formula[rows3]
group3_theta_pval = dca_pvalues$theta_formula[rows3] 

group4_mean_pval = dca_pvalues$mean_formula[rows4]
group4_theta_pval = dca_pvalues$theta_formula[rows4] 

char1 = as.character(round(deseq2_qcut, digits = 1))
char2 = as.character(round(dca_direct_qcut, digits = 1))
  

figure_filename = paste0("figures/step10k_DCA_log_formula_four_groups_pvalue_hist_", 
                         char1, "_", char2, ".pdf")
pdf(figure_filename, width=12, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,3,1))

hist(log10(group1_mean_pval), 
     main="DESeq2 sign, dca_direct sign\nmean", 
     xlab="log10(p-value)", seq(-6.5, 0, by=0.5))
hist(log10(group2_mean_pval), 
     main="DESeq2 nosign, dca_direct sign\nmean", 
     xlab="log10(p-value)", seq(-6.5, 0, by=0.5))
hist(group3_mean_pval, xlim = c(0, 1),
     main="DESeq2 sign, dca_direct nosign\nmean", 
     xlab="p-value", breaks = 20)
hist(group4_mean_pval, xlim = c(0, 1),
     main="DESeq2 nosign, dca_direct nosign\nmean", 
     xlab="p-value", breaks = 20)  

hist(log10(group1_theta_pval), 
     main="DESeq2 sign, dca_direct sign\npseudo theta", 
     xlab="log10(p-value)", seq(-6.5, 0, by=0.5))
hist(log10(group2_theta_pval), 
     main="DESeq2 nosign, dca_direct sign\npseudo theta", 
     xlab="log10(p-value)", seq(-6.5, 0, by=0.5))
hist(group3_theta_pval, xlim = c(0, 1),
     main="DESeq2 sign, dca_direct nosign\npseudo theta", 
     xlab="p-value", breaks = 20)
hist(group4_theta_pval, xlim = c(0, 1),
     main="DESeq2 nosign, dca_direct nosign\npseudo theta", 
     xlab="p-value", breaks = 20)  

dev.off()

  
  

settings = c("DESeq2 sign, dca sign", 
             "DESeq2 nosign, dca sign", 
             "DESeq2 sign, dca nosign", 
             "DESeq2 nosign, dca nosign")

mean_pvalue_list = list()

mean_pvalue_list[[1]] = group1_mean_pval
mean_pvalue_list[[2]] = group2_mean_pval
mean_pvalue_list[[3]] = group3_mean_pval
mean_pvalue_list[[4]] = group4_mean_pval


theta_pvalue_list = list()

theta_pvalue_list[[1]] = group1_theta_pval
theta_pvalue_list[[2]] = group2_theta_pval
theta_pvalue_list[[3]] = group3_theta_pval
theta_pvalue_list[[4]] = group4_theta_pval
  
get_prop_001 <- function(vec){
  proportion = round(sum(vec < 0.001)/length(vec), digits = 3)
  return(format(proportion, nsmall = 3))
}

get_prop_01 <- function(vec){
  proportion = round(sum(vec < 0.01)/length(vec), digits = 3)
  return(format(proportion, nsmall = 3))
}

props_mean_001 = sapply(mean_pvalue_list, get_prop_001)
props_theta_001  = sapply(theta_pvalue_list, get_prop_001)

props_mean_01 = sapply(mean_pvalue_list, get_prop_01)
props_theta_01  = sapply(theta_pvalue_list, get_prop_01)

ranksum_mean_mat = matrix(NA, ncol = 4, nrow = 4)
ranksum_theta_mat = matrix(NA, ncol = 4, nrow = 4)
  
for (i in 1:3){
  for (j in (i+1):4){
    temp_pvalue = wilcox.test(mean_pvalue_list[[i]],
                        mean_pvalue_list[[j]])$p.value
    if (temp_pvalue >= 0.001){
      temp_pvalue = 
        format(round(temp_pvalue, digits = 3), nsmall = 3)
    }else{
      temp_pvalue = 
        formatC(temp_pvalue, format = "e", digits = 1)
    }
    ranksum_mean_mat[i,j] = temp_pvalue
    ranksum_mean_mat[j,i] = ranksum_mean_mat[i,j]
  }
}

for (i in 1:3){
  for (j in (i+1):4){
    temp_pvalue = wilcox.test(theta_pvalue_list[[i]],
                        theta_pvalue_list[[j]])$p.value
    if (temp_pvalue >= 0.001){
      temp_pvalue = 
        format(round(temp_pvalue, digits = 3), nsmall = 3)
    }else{
      temp_pvalue = 
        formatC(temp_pvalue, format = "e", digits = 1)
    }
    ranksum_theta_mat[i,j] = temp_pvalue
    ranksum_theta_mat[j,i] = ranksum_theta_mat[i,j]
  }
} 

  
# keep the ngene part just for verification
ngenes_vec = c(length(group1_mean_pval), 
               length(group2_mean_pval), 
               length(group3_mean_pval), 
               length(group4_mean_pval))

df_four_groups = data.frame("setting" = settings, 
                            "prop_less_than0.001_mean" = props_mean_001,
                            "prop_less_than0.001_theta" = props_theta_001,
                            "prop_less_than0.01_mean" = props_mean_01,
                            "prop_less_than0.01_theta" = props_theta_01,
                            "ngenes_in_group" = ngenes_vec)

write.csv(df_four_groups, 
          file = paste0("res/step10k_log_formula_four_groups_", 
                        char1, "_", char2, 
                        "_proportion.csv"),
          row.names = FALSE)

  
#group_names = c("pp", "np", "pn", "nn")

df_ranksum_mean = cbind(settings, ranksum_mean_mat)
colnames(df_ranksum_mean) = c("names", settings)

df_ranksum_theta = cbind(settings, ranksum_theta_mat)
colnames(df_ranksum_theta) = c("names", settings)

write.csv(df_ranksum_mean, 
          file = paste0("res/step10k_log_formula_four_groups_", 
                        char1, "_", char2, 
                        "_ranksum_pvalues_mean.csv"),
          row.names = FALSE)

write.csv(df_ranksum_theta, 
          file = paste0("res/step10k_log_formula_four_groups_", 
                        char1, "_", char2, 
                        "_ranksum_pvalues_theta.csv"),
          row.names = FALSE)



# further split each group into four groups
#   -- both mean and theta pvalue < 0.001
#   -- mean pvalue < 0.001, theta pvalue >= 0.001
#   -- mean pvalue >= 0.001, theta pvalue < 0.001
#   -- mean pvalue >= 0.001, theta pvalue >= 0.001

further_four_groups_mat = matrix(NA, ncol = 4, nrow = 4)

for (i in 1:4){
  further_four_groups_mat[i, 1] = 
    sum((mean_pvalue_list[[i]]<0.001) & (theta_pvalue_list[[i]]<0.001))/length(mean_pvalue_list[[i]])
  further_four_groups_mat[i, 2] = 
    sum((mean_pvalue_list[[i]]<0.001) & (!(theta_pvalue_list[[i]]<0.001)))/length(mean_pvalue_list[[i]])
  further_four_groups_mat[i, 3] = 
    sum((!(mean_pvalue_list[[i]]<0.001)) & (theta_pvalue_list[[i]]<0.001))/length(mean_pvalue_list[[i]])
  further_four_groups_mat[i, 4] = 
    sum((!(mean_pvalue_list[[i]]<0.001)) & (!(theta_pvalue_list[[i]]<0.001)))/length(mean_pvalue_list[[i]])
  }

df_further = cbind(paste0(settings, 
                          rep(" (", 4), 
                          sapply(mean_pvalue_list, length),
                          rep(")", 4)),
                   round(further_four_groups_mat, digits = 3))

colnames(df_further)[2:5] = c("mean sign, theta sign", 
                              "mean sign, theta nosign",
                              "mean nosign, theta sign",
                              "mean nosign, theta nosign")

write.csv(df_further, 
          file = paste0("res/step10k_further_four_groups_", char1, "_", char2, ".csv"), 
          row.names = FALSE)
  
  
gc()

sessionInfo()
q(save="no")