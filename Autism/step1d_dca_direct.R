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

grp1 = gsub("PFC_", "", grp)
grp1



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

library(ideas)



# number of cores for multi-core computation
nCore = 12
#nCore = Sys.getenv("SLURM_CPUS_ON_NODE")

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")
 
# The outpuf of DCA are too large to save at GitHub, e.g., 
# -rw-rw---- 1  1383783222 Dec 21 02:00 L2_3_mean_norm.tsv
# -rwxrwx--- 1   544297480 Dec 16 22:37 L2_3_mean.tsv.gz
# -rwxrwx--- 1   486210994 Dec 16 22:38 L2_3_pi.tsv.gz
# so we access them from this local folder:
data.dca.dir = "data/dca_PFC_all/"
  
data.dir  = "data/"

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
class(dat1)

dim(dat1)
dat1[1:5,1:4]


# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

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

genes = rownames(dat1)

# ------------------------------------------------------------------------
# add read-depth information
# ------------------------------------------------------------------------

dim(meta_cell)
meta_cell[1:2,]

table(meta_cell$cell == colnames(dat1))

rd_cell = colSums(dat1)
summary(rd_cell)

meta_cell$rd = rd_cell

dim(meta_cell)
meta_cell[1:2,]

names(meta_cell)[1] = "cell_id"
mean(meta_cell$cell_id == colnames(dat1))


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
  
dca_par_list = list(dca_mean, dca_disp, dca_pi)

gc()

rm(dca_mean)
rm(dca_disp)
rm(dca_pi)
rm(dat1)
gc()
gc()
gc()

# ========================================================================
# perform testing
# ========================================================================

set.seed(2020)
  
date()
  
  
# ---------------------------------------------------------------
# estimate distance across individuals
# ---------------------------------------------------------------


var2test      = "diagnosis"
var2adjust    = c("age", "sex", "Seqbatch", "RIN")

var2test_type = "binary"
var_per_cell  = c("rd")


date()

dist_list = list()

for(fit_method in c("dca_direct")){
  for(d_metric in c("Was", "JSD")){
    message(sprintf("fit_method: %s, d_metric: %s\n", fit_method, d_metric))
    message(date())
    
    label = paste(fit_method, d_metric, sep="_")
    
    dist1 = ideas_dist(dca_par_list, meta_cell, meta_ind, 
                       var_per_cell, var2test, var2test_type,
                       d_metric = d_metric, 
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
# pval calculation by kernel regression, ZINB
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
# pval calculation by permanova
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

gh = list()
for(k in 2:ncol(df)){
  method_nm = names(df)[k]
  gh[[k-1]] = ggplot(df, aes_string(x = method_nm)) + 
    labs(title = method_nm) + 
    geom_histogram(color = "darkblue", fill = "lightblue", 
                   breaks = seq(0,1,by = 0.02))
}

fig.name  = sprintf("figures/step1d_dca_direct_pval_hist_%s.pdf", 
                    grp)
file.name = sprintf("res/step1d_dca_direct_pvals_%s.tsv", 
                    grp)


pdf(fig.name, width = 9, height = 9)
ggarrange(plotlist=gh, ncol = 2, nrow = 2)
dev.off()

fwrite(df, file = file.name, sep = "\t")

gc()

sessionInfo()
q(save="no")
