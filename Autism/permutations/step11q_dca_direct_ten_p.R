# modified based on the original step11q_dca_direct_ten_p.R
# to deal with the case when PS only has one gene to rerun

# step11q_dca_direct_ten_p.R
# the main part of the code is carried over from 
# step11i_dca_direct_p.R
# except that this time we do 10 permutations and combine the results
#   - put results from 10 permutations into list




# majority of the code is carried over from step11d_dca_direct.R
# add the permutation part
# change the step name and add "_p" to files to save


# this version is based on step10b, different from the step10b version 
# on cluster in the sense that var_per_cell is set to rd, although
# this does not make any difference
# may remove this input item in improved version
# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

grp
# grp1 is used to locate DCA results
grp1 = gsub("PFC_", "", grp)
grp1


#reRUN = FALSE

# ========================================================================
# read input
# ========================================================================

library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(doParallel)
library(doRNG)
library(svd)
#library(ideas)
library(MiRKAT)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)


library(pscl)
library(emdbook)
library(foreach)
library(stats)

#source("./ideas_2/R/ideas_dist.R")
#source("./ideas_2/R/fit_zinb.R")
#source("./ideas_2/R/fit_nb.R")
#source("./ideas_2/R/divergence.R")
#source("./ideas_2/R/permanova_utilities.R")
#source("./ideas_2/R/permanova.R")

library(ideas)

print(Sys.getenv("R_HOME"))

# number of cores for multi-core computation
nCore = Sys.getenv("SLURM_CPUS_ON_NODE")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")
 
# The outpuf of DCA are too large to save at GitHub, e.g., 
# -rw-r--r--  1 wsun  staff   549M Mar 21 21:36 L2_3_dispersion.tsv.gz
# -rw-r--r--  1 wsun  staff   519M Mar 21 21:34 L2_3_mean.tsv.gz
# -rw-r--r--  1 wsun  staff   464M Mar 21 21:35 L2_3_pi.tsv.gz
# so we access them from this local folder:
data.dir  = "./data/dca_PFC_all"

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

# change to connect to the original dat1
# compare with step4m
dat1 = readRDS(file.path(sprintf("res/step11_count_matrix_%s.rds", grp)))

class(dat1)

dim(dat1)
dat1[1:5,1:4]

genes = rownames(dat1)
genes[1:5]


meta_ind = fread(file=sprintf("res/step11_meta_ind_%s.tsv", grp), 
                 sep="\t", data.table = FALSE)

meta_cell = fread(file=sprintf("res/step11_meta_cell_%s.tsv", grp), 
                  sep="\t", data.table = FALSE)

dim(meta_ind)
dim(meta_cell)

head(meta_ind)
summary(meta_ind)
meta_cell[1:2,]


names(meta_cell)[1] = "cell_id"

mean(meta_cell$cell_id == colnames(dat1))


# ------------------------------------------------------------------------
# read in DCA estimates
# ------------------------------------------------------------------------
  
f_mean = file.path(data.dir, paste0(grp1, "_mean_norm.tsv"))
f_disp = file.path(data.dir, paste0(grp1, "_dispersion.tsv.gz"))
f_pi   = file.path(data.dir, paste0(grp1, "_pi.tsv.gz"))
  
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
  
dca_par_list = list(dca_mean, dca_disp, dca_pi)
gc()
  

  
date()
  

gc()
  
rm(dca_mean)
rm(dca_disp)
rm(dca_pi)
rm(dat1)
gc()
gc()
gc()
  
  
summary(meta_cell$mitoPercent)
summary(meta_cell$riboPercent)
  
min(meta_cell$mitoPercent[meta_cell$mitoPercent>0])
min(meta_cell$riboPercent[meta_cell$riboPercent>0])
  
table(meta_cell$mitoPercent == 0)
table(meta_cell$riboPercent == 0)

meta_cell$mitoPercent1 = meta_cell$mitoPercent + 0.005
  

  
  
var2test      = "diagnosis"
var2adjust    = c("age", "sex", "Seqbatch", "RIN")

var2test_type = "binary"
#var_per_cell  = c("rd", "riboPercent", "mitoPercent1")
var_per_cell  = c("rd")




# prepare for permuting the labels

ori_diagnosis = meta_ind$diagnosis

# prepare a list to hold the permutation results
df_list = list()

# set random seeds
seed_v = c(5608, 9903, 9968, 7527, 7879, 
           3760, 2066, 7577, 5926, 9670)


for (j in 1:10){
  
  set.seed(seed_v[j])
  print(paste("j = ", j))
  print(paste("random seed for permutation = ", seed_v[j]))
  
  meta_ind$diagnosis = sample(ori_diagnosis)
  print(meta_ind$diagnosis)
  head(meta_ind)
  summary(meta_ind)
  
  set.seed(2020) # this should not affect anything
  
  date()

  dist_list = list()

  for(fit_method in c("dca_direct")){
    for(d_metric in c("Was", "JSD")){
      message(sprintf("fit_method: %s, d_metric: %s\n", fit_method, d_metric))
      message(date())
    
      label = paste(fit_method, d_metric, sep="_")
      
      dist1 = ideas_dist(dca_par_list, meta_cell, meta_ind, 
                         var_per_cell, var2test, var2adjust, 
                         var2test_type, d_metric = d_metric, 
                         fit_method = fit_method)
      dist_list[[label]] = dist1
    }
  }

  date()

  # not sure why, it will take a few gc() to see the reduction of memory usage
  for(gi in 1:10){
    gc()
  }
  gc()

  # ---------------------------------------------------------------
  # STEP 2a: pval calculation by kernel regression, ZINB
  # ---------------------------------------------------------------

  y = as.numeric(as.factor(meta_ind$diagnosis)) - 1
  table(y)

  X = model.matrix(~ age + sex + Seqbatch + RIN, data=meta_ind)

  dim(X)
  X[1:2,]
  X = X[,-1]

  n_gene = nrow(dist_list[[1]])
  pval_KR = matrix(NA, nrow = n_gene, ncol = length(dist_list))
  rownames(pval_KR) = genes
  colnames(pval_KR) = paste("KR", names(dist_list), sep="_")


  set.seed(905)

  date()
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
                    method = "permutation", nperm = 4999)
        pval = m1$p_values
      
        if(pval < 0.1){
          m1 = MiRKAT(y = y, X = X, Ks = Ki, out_type = "D", 
                      method = "permutation", nperm = 9999)
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

  n_perm  = 4999
  n_perm2 = 9999
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
    pval_k  = permanova(dist_k, meta_ind, var2test, var2adjust, 
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
      
      pval_kr = permanova(rerun_dist_k, meta_ind, var2test, var2adjust, 
                          var2test_type, n_perm=n_perm2, r.seed=r.seed)
      
      #pval_kr = permanova(dist_k[w2rerun,,], meta_ind, var2test, var2adjust, 
      #                    var2test_type, n_perm=n_perm2, r.seed=r.seed)
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
  head(df)
  
  df_list[[j]] = df



}

file.name = sprintf("res/step11q_dca_direct_pvals_%s_ten_p.rds", 
                    grp)


saveRDS(df_list, file.name)

gc()

sessionInfo()
q(save="no")
