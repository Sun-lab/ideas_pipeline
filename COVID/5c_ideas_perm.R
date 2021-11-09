# compared to 3c, this version increases the time of permuation
# to 49999 and 99999


# compared to 1c, this version considers more genes
# by setting filtering criterion to keep genes
# appearning in at last 10% of the cells


# ideas with nb, read depth as the only covariates in variable per cell

# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'CD8+Tcells_1' as default.\n")
  grp = "CD8+Tcells_1"
}else{
  eval(parse(text=args[[1]]))
}

grp


print(Sys.getenv("R_HOME"))

# ========================================================================
# libraries and path
# ========================================================================

library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(doParallel)
library(doRNG)
library(svd)
library(MiRKAT)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_bw())


# no need to load these libraries again here if load ideas
library(pscl)
library(emdbook)
library(foreach)
library(stats)

library(ideas)



# number of cores for multi-core computation
nCore = 4
#nCore = Sys.getenv("SLURM_CPUS_ON_NODE")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")

data.dir = "../../ideas_data/COVID/PBMC_10x"
label.dir = "../../ideas_data/COVID"

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


# ------------------------------------------------------------------------
# add read-depth information
# ------------------------------------------------------------------------

rd_cell = colSums(dat1)
summary(rd_cell)

meta2kp$rd = rd_cell


# ------------------------------------------------------------------------
# read in and attach perm information
# ------------------------------------------------------------------------

perm_table = fread(file.path(label.dir, "1b_permlabel.tsv"),header=TRUE)
dim(perm_table)
perm_table[1:2,1:3]

# check that each individual has one label
table(tapply(perm_table$p1, meta2kp$individual, 
             function(v){length(unique(v))}))

# get label after permutation
iperm = 1
cur_perm_cell_label = data.frame(perm_table)[, iperm]
cur_perm_ind_label = rep(NA, nrow(meta_ind))

for (i in 1:nrow(meta_ind)){
  wi = which(meta2kp$individual == meta_ind$individual[i])
  cur_perm_ind_label[i] = cur_perm_cell_label[wi[1]]
}


meta_ind$diagnosis = cur_perm_ind_label

# only keep the 14 individuals with mild or severe label
# after permutation
ind2kp = which(cur_perm_ind_label != "")
meta_ind_perm = meta_ind[ind2kp, ]
table(meta_ind_perm$diagnosis)

# only keep the cells corresponding to these 14 individuals
cell2kp_perm = which(meta2kp$individual %in% meta_ind_perm$individual)
meta2kp_perm = meta2kp[cell2kp_perm,]
dim(meta2kp_perm)
dat1_perm = dat1[, cell2kp_perm]
dim(dat1_perm)

table(colnames(dat1_perm) == meta2kp_perm$cell_id)


# ---------------------------------------------------------------
# estimate distance across individuals
# ---------------------------------------------------------------



count_matrix_perm  = as.matrix(dat1_perm)

var2test      = "diagnosis"
var2adjust    = c("age", "sex")

var2test_type = "binary"
var_per_cell  = c("rd")


dim(count_matrix_perm)
count_matrix_perm[1:2,1:4]


genes = rownames(count_matrix_perm)
rm(dat1)
rm(dat1_perm)
gc()
gc()

set.seed(2020)

date()

dist_list = list()

for(fit_method in c("nb")){
  for(d_metric in c("Was", "JSD")){
    message(sprintf("fit_method: %s, d_metric: %s\n", fit_method, d_metric))
    message(date())
    
    label = paste(fit_method, d_metric, sep="_")
    dist1 = ideas_dist(count_matrix_perm, meta2kp_perm, meta_ind_perm, 
                         var_per_cell, var2test, var2test_type, 
                         d_metric = d_metric, 
                         fit_method = fit_method)
    dist_list[[label]] = dist1
  }
}

date()

lapply(dist_list, dim)

dist_list$nb_Was[1,1:3,1:3]
dist_list$nb_JSD[1,1:3,1:3]

rm(count_matrix_perm)
gc()
gc()

# ---------------------------------------------------------------
# STEP 2a: pval calculation by kernel regression, ZINB
# ---------------------------------------------------------------


y = as.numeric(as.factor(meta_ind_perm$diagnosis)) - 1
table(y)

X = model.matrix(~ age + sex, data=meta_ind_perm)

dim(X)
X[1:2, ]
X = X[,-1]

n_gene = nrow(dist_list[[1]])
pval_KR = matrix(NA, nrow = n_gene, ncol = length(dist_list))
rownames(pval_KR) = genes
colnames(pval_KR) = paste("KR", names(dist_list), sep="_")

n_perm  = 49999
n_perm2 = 99999
set.seed(905)

for(k in 1:length(dist_list)){
  message(names(dist_list)[k])
  message(date())
  dist_k  = dist_list[[k]]
  
  pval_KR[,k] = foreach(i_g = 1:dim(dist_k)[1], .combine = "c") %dorng% {

    Di = dist_k[i_g,,]
    
    if(any(is.na(Di))){
      pval = NA 
    }else{
      Ki = D2K(Di)
      m1 = MiRKAT(y = y, X = X, Ks = Ki, out_type = "D", 
                  method = "permutation", nperm = n_perm)
      pval = m1$p_values
      
      if(pval < 0.1){
        m1 = MiRKAT(y = y, X = X, Ks = Ki, out_type = "D", 
                    method = "permutation", nperm = n_perm2)
        pval = m1$p_values
      }
    }
    pval
  }
}

date()

dim(pval_KR)
pval_KR[1:2,]


# ---------------------------------------------------------------
# STEP 2b: pval calculation by permanova
# ---------------------------------------------------------------

n_perm  = 49999
n_perm2 = 99999
r.seed  = 904

pval_PS = matrix(NA, nrow=n_gene, ncol=length(dist_list))
rownames(pval_PS) = genes
colnames(pval_PS) = names(dist_list)
colnames(pval_PS) = paste("PS", names(dist_list), sep="_")

date()

for(k in 1:length(dist_list)){
  message(names(dist_list)[k])
  message(date())
  dist_k  = dist_list[[k]]
  pval_k  = permanova(dist_k, meta_ind_perm, var2test, var2adjust, 
                      var2test_type, n_perm=n_perm, r.seed=r.seed)
  w2rerun = which(pval_k < 0.1)
  if(length(w2rerun) > 0){
    if (length(w2rerun) == 1){
      slice_dist_k = dist_k[w2rerun, , ]  
      rerun_dist_k = array(dim = c(1, dim(slice_dist_k)))
      rerun_dist_k[1, , ] = slice_dist_k
    }else{
      rerun_dist_k = dist_k[w2rerun, , ]
    }
    pval_kr = permanova(rerun_dist_k, meta_ind_perm, var2test, var2adjust, 
                        var2test_type, n_perm=n_perm2, r.seed=r.seed)
    pval_k[w2rerun] = pval_kr
  }
  pval_PS[,k] = pval_k
}

date()

summary(pval_KR)
summary(pval_PS)

# ------------------------------------------------------------------------
# summarize and save the results
# ------------------------------------------------------------------------


df = data.frame(gene=genes, pval_KR, pval_PS)
dim(df)
df[1:5,]

gh = list()
for(k in 2:ncol(df)){
  method_nm = names(df)[k]
  gh[[k-1]] = ggplot(df, aes_string(x = method_nm)) + 
    labs(title = method_nm) + 
    geom_histogram(color = "darkblue", fill = "lightblue", 
                   breaks = seq(0,1,by = 0.02))
}


fig.name  = sprintf("figures/5c_nb_pval_hist_%s_perm.pdf", grp)
file.name = sprintf("res/5c_nb_pvals_%s_perm.tsv", grp)


pdf(fig.name, width = 10, height = 10)
ggarrange(plotlist=gh, ncol = 2, nrow = 2)
dev.off()

fwrite(df, file = file.name, sep = "\t")

gc()

sessionInfo()
q(save="no")
