# modifed from step3a_DCA_formula_helper.R
# on Autism data


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
data.dir = "../../ideas_data/COVID/PBMC_10x"
data.dca.dir = "../../ideas_data/COVID/PBMC_10x/dca_zinb_celltypes"

grp =  "CD8+Tcells_1"



## data processing

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"), 
                  stringsAsFactors=TRUE)
dim(cell_info)
cell_info[1:2,]

# ------------------------------------------------------------------------
# read in count data of one celltype
# ------------------------------------------------------------------------

dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat)
class(dat)
dat[1:5,1:4]

# ------------------------------------------------------------------------
# read in covid donor information
# ------------------------------------------------------------------------

covid_donor_info = 
  read.csv(file.path(data.dir, "covid_donor_info_from_mmc1.csv"), 
           header = TRUE)

dim(covid_donor_info)
covid_donor_info[1:2,]
summary(covid_donor_info)

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat) %in% cell_info$cell)

meta = cell_info[match(colnames(dat), cell_info$cell),]
dim(meta)
meta[1:2,]

summary(meta)

meta$donor = as.factor(meta$donor)

summary(meta$nCount_RNA/meta$nFeature_RNA)


# filter out cells from control samples
meta_covid = meta[which(meta$group_per_sample != "control"),]
dim(meta_covid)

table(meta_covid$group_per_sample)
table(meta_covid$disease_stage)
table(meta_covid$donor)

table(meta_covid$donor, meta_covid$group_per_sample)

df_donor = as.data.frame(table(meta_covid$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 10)]

meta2kp = meta_covid[which(meta_covid$donor %in% donor2kp),]
dim(meta2kp)
meta2kp[1:2,]
table(meta2kp$donor)
length(unique(meta2kp$donor))

cell2kp_index = which(meta$cell %in% meta2kp$cell)

# select counts in the cells to keep
dat1 = dat[, cell2kp_index]
table(colnames(dat1) == meta2kp$cell)


# check how many samples (each sample contains multiple cells) each donor has
table(tapply(meta2kp$sampleID, meta2kp$donor, function(v){length(unique(v))}))

sort(table(paste(meta2kp$donor, meta2kp$group_per_sample, sep=":")))

# adjust certain column names in meta2kp to match the requirement 
# of ideas and for the ease of later processing
colnames_meta2kp = names(meta2kp)
names(meta2kp)[which(colnames_meta2kp=="cell")] = "cell_id"
names(meta2kp)[which(colnames_meta2kp=="donor")] = "individual"

names(meta2kp)[which(colnames_meta2kp=="group_per_sample")] = "diagnosis"

# ------------------------------------------------------------------------
# generate individual level information
# ------------------------------------------------------------------------

meta2kp$diagnosis = droplevels(meta2kp$diagnosis)
table(meta2kp$diagnosis)

meta2kp$individual = droplevels(meta2kp$individual)
meta2kp$sex = droplevels(meta2kp$sex)

meta_ind = distinct(meta2kp[,c('individual', 'diagnosis', 'sex')])
table(meta_ind$diagnosis, meta_ind$sex)

# add exact age information
donor_info_match = 
  covid_donor_info[match(meta_ind$individual, covid_donor_info$donor),]
# double check that the condition and sex features match
table(donor_info_match$condition == meta_ind$diagnosis)
table(which(donor_info_match$sex=="f") == which(meta_ind$sex=="female"))
meta_ind$age = donor_info_match$age

sort(meta_ind$age[which(meta_ind$diagnosis=="mild")])
sort(meta_ind$age[which(meta_ind$diagnosis=="severe")])

meta_ind$age = scale(meta_ind$age)

table(meta_ind$diagnosis)


if(nrow(meta_ind) != length(unique(meta2kp$individual))){
  stop("there is non-unique information\n")
}

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

w2kp = which(n.zeros < 0.9*ncol(dat1))
dat1 = dat1[w2kp,]

dim(dat1)
dat1[1:5,1:4]


genes = rownames(dat1)



# ------------------------------------------------------------------------
# add read-depth information
# ------------------------------------------------------------------------

rd_cell = colSums(dat1)
summary(rd_cell)

meta2kp$rd = rd_cell



# ------------------------------------------------------------------------
# read in DCA estimates
# ------------------------------------------------------------------------

f_mean = file.path(data.dca.dir, paste0(grp, "_mean_norm.tsv"))
f_disp = file.path(data.dca.dir, paste0(grp, "_dispersion.tsv"))
f_pi   = file.path(data.dca.dir, paste0(grp, "_pi.tsv"))

dca_mean = fread(f_mean, sep="\t", data.table = FALSE)
dca_disp = fread(f_disp, sep="\t", data.table = FALSE)
dca_pi   = fread(f_pi,sep="\t", data.table = FALSE)

dim(dca_mean)
dim(dca_disp)
dim(dca_pi)

dca_mean[1:2,1:5]
dca_disp[1:2,1:5]
dca_pi[1:2,1:5]

table(meta$cell == colnames(dca_mean)[-1])
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
dca_disp = data.matrix(dca_disp[,-1])
dca_pi   = data.matrix(dca_pi[,-1])

dca_mean = dca_mean[, cell2kp_index]
dca_disp = dca_disp[, cell2kp_index]
dca_pi   = dca_pi[, cell2kp_index]

table(colnames(dat1) == colnames(dca_mean))
table(colnames(dat1) == colnames(dca_disp))
table(colnames(dat1) == colnames(dca_pi))


gc()
gc()
gc()
gc()



n_gene   = nrow(dca_mean)

mean_matrix = matrix(nrow = n_gene, ncol = length(meta_ind$individual))
var_matrix = matrix(nrow = n_gene, ncol = length(meta_ind$individual))


for (i_g in 1:n_gene){
  
  mean_v = rep(NA, nrow(meta_ind))
  var_v = rep(NA, nrow(meta_ind))
  
  for (j in 1:nrow(meta_ind)){
    cur_ind = meta_ind$individual[j]  
    cell_columns = which(meta2kp$individual == cur_ind)
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
  
saveRDS(mean_matrix, file = "res/7a_DCA_formula_mean_matrix.rds")
saveRDS(var_matrix, file = "res/7a_DCA_formula_var_matrix.rds")





gc()

sessionInfo()
q(save="no")