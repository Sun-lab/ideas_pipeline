# step11p_ideas_ten_p.R
# the main of the code is carried over from step11h_ideas_p.R
#    -  modify the output data frame to hold results from 
#       ten permutations
#    -  put results from ten permutations into list





# the main part of code is carried over from step11c_ideas.R
# add the permutation part
# change the step name and add "_p" to files to save


# formalize ideas with nb only on raw data
# fixed random seed
# with RNGkind("L'Ecuyer-CMRG")
# based on:
# step4b_ideas_m_nb.R
# step10c_ideas_use_DCA_nb_seed_test.R, which is based on
# step9r_ideas_use_DCA_nb_seed_m.R

# this version uses rd only
# but the rd only notation in the file names are omitted
# to formalize the files
# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
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
#library(ideas)
library(MiRKAT)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_bw())





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

# number of cores for multi-core computation
nCore = Sys.getenv("SLURM_CPUS_ON_NODE")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

dat1 = readRDS(file.path(sprintf("res/step11_count_matrix_%s.rds", grp)))
class(dat1)

dim(dat1)
dat1[1:5,1:4]

# ------------------------------------------------------------------------
# read cell level and individual level information
# ------------------------------------------------------------------------

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
table(meta_cell$cell_id == colnames(dat1))





count_matrix  = as.matrix(dat1)

var2test      = "diagnosis"
var2adjust    = c("age", "sex", "Seqbatch", "RIN")


var2test_type = "binary"
#var_per_cell  = c("rd", "riboPercent", "mitoPercent1")
var_per_cell  = c("rd")

summary(meta_cell$mitoPercent)
summary(meta_cell$riboPercent)

min(meta_cell$mitoPercent[meta_cell$mitoPercent>0])
min(meta_cell$riboPercent[meta_cell$riboPercent>0])

table(meta_cell$mitoPercent == 0)
table(meta_cell$riboPercent == 0)

meta_cell$mitoPercent1 = meta_cell$mitoPercent + 0.005

dim(count_matrix)
count_matrix[1:2,1:4]



genes = rownames(count_matrix)
rm(dat1)
gc()
gc()



# prepare for permuting the labels

ori_diagnosis = meta_ind$diagnosis

# prepare a list to hold the permutation results
df_list = list()

# set random seeds
seed_v = c(5608, 9903, 9968, 7527, 7879, 
           3760, 2066, 7577, 5926, 9670)

# loop 10 times
for (j in 1:10){
  
  set.seed(seed_v[j])
  print(paste("j = ", j))
  print(paste("random seed for permutation = ", seed_v[j]))
  
  meta_ind$diagnosis = sample(ori_diagnosis)
  print(meta_ind$diagnosis)
  head(meta_ind)
  summary(meta_ind)



# ---------------------------------------------------------------
# estimate distance across individuals
# ---------------------------------------------------------------


  set.seed(2020)

  date()

  dist_list = list()

  for(fit_method in c("nb")){
    for(d_metric in c("Was", "JSD")){
      message(sprintf("fit_method: %s, d_metric: %s\n", fit_method, d_metric))
      message(date())
    
      label = paste(fit_method, d_metric, sep="_")
      dist1 = ideas_dist(count_matrix, meta_cell, meta_ind, 
                         var_per_cell, var2test, var2adjust, 
                         var2test_type, d_metric = d_metric, 
                         fit_method = fit_method)
      dist_list[[label]] = dist1
      # comment out the saving line below to avoid saving too many files
      # saveRDS(dist1, fnm) 
    }
  }

  date()

  lapply(dist_list, dim)

  dist_list$nb_Was[1,1:3,1:3]
  dist_list$nb_JSD[1,1:3,1:3]

  # ---------------------------------------------------------------
  # STEP 2a: pval calculation by kernel regression, ZINB
  # ---------------------------------------------------------------


  y = as.numeric(as.factor(meta_ind$diagnosis)) - 1
  table(y)

  X = model.matrix(~ age + sex + Seqbatch + RIN, data=meta_ind)


  dim(X)
  X[1:2, 1:5]
  X = X[,-1]

  n_gene = nrow(dist_list[[1]])
  pval_KR = matrix(NA, nrow = n_gene, ncol = length(dist_list))
  rownames(pval_KR) = genes
  colnames(pval_KR) = paste("KR", names(dist_list), sep="_")

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

    dim(dist_k)
    pval_k  = permanova(dist_k, meta_ind, var2test, var2adjust, 
                        var2test_type, n_perm=n_perm, r.seed=r.seed)
  
    w2rerun = which(pval_k < 0.1)
    if(length(w2rerun) > 0){
      pval_kr = permanova(dist_k[w2rerun,,], meta_ind, var2test, var2adjust, 
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

  df_list[[j]] = df

}


rm(count_matrix)
gc()
gc()


file.name = sprintf("res/step11p_nb_pvals_%s_ten_p.rds", grp)

saveRDS(df_list, file.name)


gc()

sessionInfo()
q(save="no")
