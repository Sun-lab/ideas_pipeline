# This code does log(mean) and log(variance) test for all 8260 genes


## current version handles L2_3 only




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
data.dca.dir = "../../ideas_data/Autism/dca_PFC_all"

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


full_dat1 = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
full_genes = rownames(full_dat1)  



dim(full_dat1)
full_dat1[1:5,1:4]

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

dat1 = full_dat1

table(colnames(dat1) %in% cell_info$cell)

meta_cell = cell_info[match(colnames(dat1), cell_info$cell),]
dim(meta_cell)
meta_cell[1:2,]

names(meta_cell)[11:12] = c("PMI", "RIN")
names(meta_cell)[15:16] = c("mitoPercent", "riboPercent")
dim(meta_cell)
meta_cell[1:2,]

summary(meta_cell)
meta_cell$age = scale(meta_cell$age)
meta_cell$PMI = scale(meta_cell$PMI)
meta_cell$RIN = scale(meta_cell$RIN)
meta_cell$Capbatch = droplevels(meta_cell$Capbatch)
meta_cell$Seqbatch = droplevels(meta_cell$Seqbatch)
meta_cell$individual = as.factor(meta_cell$individual)
summary(meta_cell)

table(meta_cell$Capbatch, meta_cell$Seqbatch)
summary(meta_cell$UMIs/meta_cell$genes)

# make sure each individual has a unique sample
table(tapply(meta_cell$sample, meta_cell$individual, 
             function(v){length(unique(v))}))

# check each individual has a unique Capbatch
table(tapply(meta_cell$Capbatch, meta_cell$individual, 
             function(v){length(unique(v))}))

table(meta_cell$cluster)
table(meta_cell$region)
table(meta_cell$diagnosis)
sort(table(paste(meta_cell$individual, meta_cell$diagnosis, sep=":")))
tt1 = table(meta_cell$individual)
sort(tt1)

# ------------------------------------------------------------------------
# generate individual level information
# ------------------------------------------------------------------------

length(unique(meta_cell$individual))

meta_ind = distinct(meta_cell[,3:12])
dim(meta_ind)
meta_ind[1:2,]
meta_ind$diagnosis = relevel(meta_ind$diagnosis, ref="Control")
table(meta_ind$diagnosis)

if(nrow(meta_ind) != length(unique(meta_cell$individual))){
  stop("there is non-unique information\n")
}

table(meta_ind$Seqbatch, meta_ind$Capbatch)



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
# load means and variances from formula based approach
# -------------------------------------------------------------------

mean_formula = readRDS("res/step3a_DCA_formula_mean_matrix.rds")
var_formula = readRDS("res/step3a_DCA_formula_var_matrix.rds")


# ------------------------------------------------------------
# pvalues on log(mean) and log(variance) from dca formula 
# ------------------------------------------------------------

v_mean_formula_pvals = rep(NA, nrow(mean_formula))

for(i_g in 1:nrow(mean_formula)){
    
  v_means = log(mean_formula[i_g, ])
  temp_mean_df = cbind(meta_ind, v_means)
  colnames(temp_mean_df)[ncol(temp_mean_df)] = "means"
  temp_lm <- lm(means ~ age + sex + Seqbatch + RIN + diagnosis, data = temp_mean_df)
  temp_pvalue  = summary(temp_lm)$coefficients[6, 4]                     
  v_mean_formula_pvals[i_g] = temp_pvalue
    
}

 
v_var_formula_pvals = rep(NA, nrow(var_formula)) 

for(i_g in 1:nrow(var_formula)){
  
  v_vars = log(var_formula[i_g, ])
  temp_var_df = cbind(meta_ind, v_vars)
  colnames(temp_var_df)[ncol(temp_var_df)] = "vars"
  temp_lm <- lm(vars ~ age + sex + Seqbatch + RIN + diagnosis, data = temp_var_df)
  temp_pvalue  = summary(temp_lm)$coefficients[6, 4]                       
  v_var_formula_pvals[i_g] = temp_pvalue
  
}



# save pvalue results out

pvalues_combined = data.frame(gene = full_genes[w2kp], 
                              mean_formula = v_mean_formula_pvals, 
                              var_formula = v_var_formula_pvals)

write.csv(pvalues_combined, 
          file = "res/step3b_formula_covariates_pvals_all.csv",
          row.names = FALSE)





gc()

sessionInfo()
q(save="no")