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
data.dca.dir = "../../ideas_data/COVID/PBMC_10x/dca_zinb_celltypes"
data.dir = "../../ideas_data/COVID/PBMC_10x"

grp = "CD8+Tcells_1"



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

full_genes = row.names(dat)

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



# -------------------------------------------------------------------
# load means and variances from formula based approach
# -------------------------------------------------------------------

mean_formula = readRDS("res/7a_DCA_formula_mean_matrix.rds")
var_formula = readRDS("res/7a_DCA_formula_var_matrix.rds")

# ------------------------------------------------------------
# write out matrix for theta
# ------------------------------------------------------------

summary(c(var_formula) - c(mean_formula))

theta_formula = (mean_formula)^2/(var_formula - mean_formula)
dim(theta_formula)
theta_formula[1:6, 1:2]

summary(c(theta_formula))

saveRDS(theta_formula, 
        file = "res/7c_DCA_formula_pseudo_theta_matrix.rds")


# ------------------------------------------------------------
# pvalues on log(theta) from dca formula 
# ------------------------------------------------------------




v_theta_formula_pvals = rep(NA, nrow(theta_formula))

for(i_g in 1:nrow(theta_formula)){
    
  v_thetas = log(theta_formula[i_g, ])
  temp_theta_df = cbind(meta_ind, v_thetas)
  colnames(temp_theta_df)[ncol(temp_theta_df)] = "thetas"
  temp_lm <- lm(thetas ~ age + sex + diagnosis,  data = temp_theta_df)
  temp_pvalue  = summary(temp_lm)$coefficients[4, 4]                     
  v_theta_formula_pvals[i_g] = temp_pvalue
    
}

 
summary(v_theta_formula_pvals)

# save pvalue results out

pvalues_combined = data.frame(gene = full_genes[w2kp], 
                              theta_formula = v_theta_formula_pvals)

write.csv(pvalues_combined, 
          file = "res/7c_formula_theta_covariates_pvals_all.csv",
          row.names = FALSE)





gc()

sessionInfo()
q(save="no")