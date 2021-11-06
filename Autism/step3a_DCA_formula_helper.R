

# This file is modified from 
# step3a_DCA_perspective.R
# to save a matrix of means and a matrix of variance 
# computed based on formulas

## current version handles L2_3 only




library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
#library(doParallel)
#library(doRNG)
#library(svd)
#library(ideas)
#library(MiRKAT)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)

library(grid)
library(gridExtra)
#big_font <- theme_grey(base_size =  16)

theme_set(theme_classic())





#setwd("~/Documents/Fred_Hutch/core_code")
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



# ------------------------------------------------------------------------
# add read-depth information
# ------------------------------------------------------------------------

dim(meta_cell)
meta_cell[1:2,]

table(meta_cell$cell == colnames(dat1))

rd_cell = colSums(dat1)
summary(rd_cell)

meta_cell$rd = rd_cell

names(meta_cell)[1] = "cell_id"
table(meta_cell$cell_id == colnames(dat1))


# ---------------------------------------------------------------------


# ----------------------------------------------------------------------------
# (2) individual plots under DCA_direct
# ----------------------------------------------------------------------------
#    - DCA_direct code

# ------------------------------------------------------------------------
# read in DCA estimates
# ------------------------------------------------------------------------

f_mean = file.path(data.dca.dir, paste0(grp1, "_mean_norm.tsv"))
f_disp = file.path(data.dca.dir, paste0(grp1, "_dispersion.tsv.gz"))
f_pi   = file.path(data.dca.dir, paste0(grp1, "_pi.tsv.gz"))

dca_mean = fread(f_mean, sep="\t", data.table = FALSE)
dca_disp = fread(f_disp)
dca_pi   = fread(f_pi)

dim(dca_mean)
dim(dca_disp)
dim(dca_pi)

dca_mean[1:2,1:5]
dca_disp[1:2,1:5]
dca_pi[1:2,1:5]

table(meta_cell$cell_id == colnames(dca_mean)[-1])
table(rownames(dat1) %in% dca_mean$V1)

w2kp = match(rownames(dat1), dca_mean$V1)

dca_mean = dca_mean[w2kp,]
dca_disp = dca_disp[w2kp,]
dca_pi   = dca_pi[w2kp,]

table(rownames(dat1) == dca_mean$V1)
table(rownames(dat1) == dca_disp$V1)
table(rownames(dat1) == dca_pi$V1)

rownames(dca_mean) = dca_mean$V1
rownames(dca_disp) = dca_disp$V1
rownames(dca_pi)   = dca_pi$V1

dca_mean = data.matrix(dca_mean[,-1])
dca_disp = data.matrix(dca_disp[,-1, with=FALSE])
dca_pi   = data.matrix(dca_pi[,-1,   with=FALSE])

table(colnames(dat1) == colnames(dca_mean))
table(colnames(dat1) == colnames(dca_disp))
table(colnames(dat1) == colnames(dca_pi))

gc()
gc()
gc()
gc()



gene_ids = rownames(dca_mean)
n_gene   = nrow(dca_mean)

mean_matrix = matrix(nrow = n_gene, ncol = length(meta_ind$individual))
var_matrix = matrix(nrow = n_gene, ncol = length(meta_ind$individual))


for (i_g in 1:n_gene){
  
  mean_v = rep(NA, nrow(meta_ind))
  var_v = rep(NA, nrow(meta_ind))
  
  for (j in 1:nrow(meta_ind)){
    cur_ind = meta_ind$individual[j]  
    cell_columns = which(meta_cell$individual == cur_ind)
    means_for_ind = c()
    vars_for_ind = c()
    for (j_c in cell_columns){
      temp_mu = as.numeric(dca_mean[i_g, j_c])
      temp_k = as.numeric(dca_disp[i_g, j_c])
      temp_pi = as.numeric(dca_pi[i_g, j_c])
      temp_mean = (1-temp_pi)*temp_mu
      # variance formula verified through simulation
      temp_var = (1-temp_pi)*(temp_mu + temp_mu^2 * (1+1/temp_k)) - (1-temp_pi)^2 * temp_mu^2
      means_for_ind = c(means_for_ind, temp_mean)
      vars_for_ind = c(vars_for_ind, temp_var)
    }
    mean_v[j] = mean(means_for_ind)
    var_v[j] = mean(vars_for_ind) + mean((means_for_ind - mean_v[j])^2)
  }
  
  mean_matrix[i_g, ] = mean_v
  var_matrix[i_g, ] = var_v
  
}
  
saveRDS(mean_matrix, file = "res/step3a_DCA_formula_mean_matrix.rds")
saveRDS(var_matrix, file = "res/step3a_DCA_formula_var_matrix.rds")





gc()

sessionInfo()
q(save="no")