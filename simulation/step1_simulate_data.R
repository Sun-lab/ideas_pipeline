#########################################################################
#                                                                       #
#                                                                       #
#                 Part I: scRNA seq data Simulation                     #
#                                                                       #
#                                                                       #
#########################################################################

# We simulate three groups of genes
#    The first group has mean expression difference,
#    The second group has variance difference.
#    The third group are equivalently expressed

# we simulate scRNAseq data per gene per cell from a zero inflated 
# negative binomial distribution. In the following codes, we 
# simulate our data based on the dataset from Velmeshev et al, 2019

# This dataset was generated using 10x Genomics platform. The read count
# data were downloaded from the link of "Gene / cell matrix (raw)" .
# from the interactive web browser at Velmeshev et al, 2019
# (https://cells.ucsc.edu/autism/).

# We apply the Deep count autoencoder(DCA) for denoising
# (https://github.com/theislab/dca) the data. Then we use the  
# output of DCA as reference to simulate our data.

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) < 5) {
  message("no enough arguments, using default values")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
  ncase    = 13      # case individuals
  nctrl    = 10      # control individuals
  ncell    = 360    # numbers of cells collected from each individuals.
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
# additional parameters
# ---------------------------------------------------------------

nCore      = 8       # number of cores for multi-core computation
nGeneMean  = 1000    # number of genes with different means in cases
nGeneVar   = 1000    # number of genes with different variance in cases
nGeneBlank = 6000    # number of genes equivalently expressed
nGeneTotal = nGeneMean + nGeneVar + nGeneBlank # total numbers of genes
nall       = ncase + nctrl

# we use the cells from one cell type (specifieid by grp1) for simulation
grp = "PFC_L2_3"
grp1 = gsub("PFC_", "", grp)
grp1

data.dir.github = "../Autism/data/"

# The outpuf of DCA are too large to save at GitHub, e.g., 
# -rw-r--r--  1 wsun  staff   549M Mar 21 21:36 L2_3_dispersion.tsv.gz
# -rw-r--r--  1 wsun  staff   519M Mar 21 21:34 L2_3_mean.tsv.gz
# -rw-r--r--  1 wsun  staff   464M Mar 21 21:35 L2_3_pi.tsv.gz
# so we access them from this local folder:
data.dir.local  = "~/research/scRNAseq/data/autism/dca_PFC_all/"

# ---------------------------------------------------------------
# initial setup
# ---------------------------------------------------------------

library(MASS)
library(Matrix)
library(emdbook)
library(moments)
library(doParallel)
library(foreach)
library(data.table)
library(pryr)
library(ggplot2)
library(ggpubr)
library(doRNG)
theme_set(theme_bw())

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

## NOTE: the data matrix will be permuted later
## so these indexes are the index in thepermuted data. 

# setting the index for the genes in three categories
i_mean  = 1:nGeneMean
i_var   = (nGeneMean + 1):(nGeneMean + nGeneVar)

set.seed(1999)

# sample gene index for genes differential expressed by mean or variance.
special_index = sample.int(nGeneTotal, (nGeneMean + nGeneVar))
mean_index    = as.numeric(special_index[i_mean])
var_index     = as.numeric(special_index[i_var])
EE_index      = (1:nGeneTotal)[-special_index]

geneIndex = list(mean_index=mean_index, var_index=var_index, 
                 EE_index=EE_index)

# ------------------------------------------------------------------------
# load real data as reference for simulation
# ------------------------------------------------------------------------

dat1 = readRDS(file.path(data.dir.github, sprintf("ct_mtx/%s.rds", grp)))
class(dat1)

dim(dat1)
dat1[1:5,1:4]

n.zeros = rowSums(dat1 == 0)
summary(n.zeros)

table(n.zeros < 0.6*ncol(dat1))
table(n.zeros < 0.8*ncol(dat1))

# only keep the top nGeneTotal genes for our simulations
w2kp = which(rank(n.zeros) <= nGeneTotal)
summary(n.zeros[w2kp]/ncol(dat1))

# ------------------------------------------------------------------------
# read in DCA output for this dataset
# ------------------------------------------------------------------------

# The major output of DCA includes 3 files, which describes the 
# mean/dipersion/dropout probabilities of each gene in each cell.
# According to their official website(https://github.com/theislab/dca):
# 
# mean.tsv represents the mean parameter of the ZINB distribution
# for each cell and gene.
# dispersion.tsv, the dispersion for each cell and gene.
# pi.tsv represent dropout probabilities for each cell and gene.

f_mean = file.path(data.dir.local, paste0(grp1, "_mean.tsv.gz"))
f_disp = file.path(data.dir.local, paste0(grp1, "_dispersion.tsv.gz"))
f_pi   = file.path(data.dir.local, paste0(grp1, "_pi.tsv.gz"))

dca_mean = fread(f_mean)
dca_disp = fread(f_disp)
dca_pi   = fread(f_pi)

dim(dca_mean)
dim(dca_disp)
dim(dca_pi)

dca_mean[1:2,1:5]
dca_disp[1:2,1:5]
dca_pi[1:2,1:5]

table(dca_mean$V1 == rownames(dat1))

t_mean = data.matrix(dca_mean[,-1, with=FALSE])
t_disp = data.matrix(dca_disp[,-1, with=FALSE])
t_drop = data.matrix(dca_pi[,-1, with=FALSE])

rownames(t_mean) = dca_mean$V1
rownames(t_disp) = dca_disp$V1
rownames(t_drop) = dca_pi$V1

t_mean = t_mean[w2kp,]
t_disp = t_disp[w2kp,]
t_drop = t_drop[w2kp,]

dim(t_mean)
t_mean[1:2,1:5]

dim(t_disp)
t_disp[1:2,1:5]

dim(t_drop)
t_drop[1:2,1:5]

summary(apply(t_drop,1,median))

gc()
rm(dca_mean)
rm(dca_disp)
rm(dca_pi)

gc()
gc()

# ------------------------------------------------------------------------
# summarize these parameters at sample level
# we want to estimate the sample_log_mean, sample_log_disp, and 
# sample_logit_drop per gene and per sample, by averaging across cells. 
# 
# In addition, we want to estimate the sd of sample_log_mean, after 
# removing the variation due to read-depth difference.
# ------------------------------------------------------------------------

col_info   = strsplit(colnames(t_mean), split="_")
sample_ids = sapply(col_info, function(x){paste(x[-1], collapse="_")})
sort(table(sample_ids))
median(sort(table(sample_ids)))

cell_rd = colSums(dat1)
cell_mean_sum = colSums(t_mean)

summary(cell_rd)
summary(cell_mean_sum)
cor(cell_rd, cell_mean_sum)

log_t_mean   = log(t(t(t_mean)*(10000/cell_mean_sum)))
logit_t_drop = log(t_drop/(1 - t_drop))

xmt = model.matrix(~log(cell_mean_sum))
dim(xmt)
xmt[1:2,]

log_t_mean_resid = matrix(NA, nrow=nrow(t_mean), ncol=ncol(t_mean))

coef = matrix(NA, nrow=nrow(log_t_mean), ncol=2)
for(i in 1:nrow(log_t_mean)){
  yi = log_t_mean[i,]
  li = lm.fit(x=xmt, y=yi)
  coef[i,] = li$coefficients
  log_t_mean_resid[i,] = li$residuals
}
summary(coef)

dim(log_t_mean_resid)
log_t_mean_resid[1:2,1:2]

tapply_median <- function(x){tapply(x, sample_ids, median)}
tapply_sd <- function(x){tapply(x, sample_ids, sd)}

sample_log_mean   = t(apply(log_t_mean,   1, tapply_median))
sample_log_disp   = t(apply(log(t_disp),  1, tapply_median))
sample_logit_drop = t(apply(logit_t_drop, 1, tapply_median))

dim(sample_log_mean)
sample_log_mean[1:2,1:3]
summary(c(sample_log_mean))

dim(sample_log_disp)
sample_log_disp[1:2,1:3]
summary(c(sample_log_disp))

dim(sample_logit_drop)
sample_logit_drop[1:2,1:3]
summary(c(sample_logit_drop))

sample_log_mean_sd = t(apply(log_t_mean_resid, 1, tapply_sd))
dim(sample_log_mean_sd)
sample_log_mean_sd[1:2,1:3]
summary(c(sample_log_mean_sd))

# the sd of log_mean across individuals: it is calculted by 
# taking the median of log_mean within individual, and then 
# calculate sd across individuals
summary(apply(sample_log_mean, 1, sd))

# the median of the sd of log_mean within individuals
summary(apply(sample_log_mean_sd, 1, median))

gc()
rm(t_mean)
rm(t_disp)
rm(t_drop)

rm(log_t_mean)
rm(log_t_mean_resid)
rm(logit_t_drop)

rm(dat1)

for(i in 1:10){
  gc()
}
gc()

# ------------------------------------------------------------------------
# check the effect of covariates
# ------------------------------------------------------------------------

meta0 = fread(file.path(data.dir.github, "meta.tsv"))
dim(meta0)
meta0[1:2,]

w2kp  = which(meta0$region == "PFC" & meta0$cluster == "L2/3")
meta0 = meta0[w2kp,]
dim(meta0)
meta0[1:2,]

meta0_ind = base::unique(meta0[,3:12])
dim(meta0_ind)
meta0_ind[1:2,]
names(meta0_ind)[9:10] = c("PMI", "RIN")

table(meta0_ind$sex, meta0_ind$diagnosis)
table(meta0_ind$Seqbatch, meta0_ind$diagnosis)

table(meta0_ind$sample == colnames(sample_log_mean))

pvals = matrix(NA, nrow=nGeneTotal, ncol=4)
colnames(pvals) = c("age", "sex", "seqBatch", "RIN")

for(i in 1:nGeneTotal){
  yi = sample_log_mean[i,]
  lmi = lm(yi ~ age, data=meta0_ind)
  pvals[i,1] = as.numeric(summary(lmi)$coefficients[2,4])
  
  lmi = lm(yi ~ sex, data=meta0_ind)
  pvals[i,2] = as.numeric(summary(lmi)$coefficients[2,4])
  
  lmi = lm(yi ~ Seqbatch, data=meta0_ind)
  pvals[i,3] = as.numeric(summary(lmi)$coefficients[2,4])
  
  lmi = lm(yi ~ RIN, data=meta0_ind)
  pvals[i,4] = as.numeric(summary(lmi)$coefficients[2,4])
}

summary(pvals)

# ------------------------------------------------------------------------
# get the number of cells per individual
# ------------------------------------------------------------------------

sort(table(paste(meta0$individual, meta0$diagnosis, sep=":")))

ncell_case = table(meta0$individual[which(meta0$diagnosis == "ASD")])
ncell_ctrl = table(meta0$individual[which(meta0$diagnosis == "Control")])
sort(ncell_case)
sort(ncell_ctrl)

# ------------------------------------------------------------------------
# simulation the 4 parameters for each gene across samples based on 
# a multivariate-normal distribution estimation for 
# log_mean, log_disp, logit_drop, log of log_mean_sd.
# ------------------------------------------------------------------------

par.names = c("mean", "dispersion", "dropout", "mean_sd")
sample_ctrl = array(dim = c(nGeneTotal, nall, 4), 
                    dimnames = list(paste0("gene", 1:nGeneTotal), 
                                    paste0("ind", 1:nall), par.names))

# first apply a normal quantile transformation to RIN, so that 
# later we can simply simulate RIN from standard normal distribution
normscore = function(vec) {
  len  = length(na.omit(vec))+1
  rank = rank(na.omit(vec))
  ties = (rank - floor(rank)) > 0
  new.vec = vec[!is.na(vec)] 
  new.vec[!ties]=qnorm(rank[!ties]/len)
  new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
  vec[!is.na(vec)] = new.vec
  vec
}

RIN.qn = normscore(meta0_ind$RIN)
RIN.simu = rnorm(nall)

set.seed(904)

for (ig in 1:nGeneTotal) {
  sample_data = cbind(
    c(sample_log_mean[ig, ]),
    c(sample_log_disp[ig, ]),
    c(sample_logit_drop[ig, ]),
    c(log(sample_log_mean_sd[ig, ]))
  )
  
  sample_data_mean = apply(sample_data, 2, mean, na.rm = TRUE)
  cov_matrix = cov(sample_data)
  
  log_mean_ig = sample_data[,1]
  lmi  = lm(log_mean_ig ~ RIN.qn)
  beta = lmi$coefficients
  
  # add some extra variance for the mean parameter
  e1 = rnorm(nall, mean=0, sd=sqrt(cov_matrix[1,1]))
  log_mean_ig_simu = beta[1] + beta[2]*RIN.simu + e1
  
  for(j in 1:nall){
    sample_data_mean_j    = sample_data_mean
    sample_data_mean_j[1] = log_mean_ig_simu[j]
    sample_ctrl[ig, j, ]  = exp(mvrnorm(1, mu = sample_data_mean_j, 
                                      Sigma = cov_matrix))  
  }
}

# double check the p-value of the covariate RIN
pvals.check = rep(NA, nGeneTotal)
for (ig in 1:nGeneTotal) {
  y_ig = log(sample_ctrl[ig, , 1])
  lm2  = lm(y_ig ~ RIN.simu)
  pvals.check[ig] = summary(lm2)$coef[2,4]
}
summary(pvals.check)

# the dropout
sample_ctrl[, , 3] = sample_ctrl[, , 3] / (1 + sample_ctrl[, , 3])

dim(sample_ctrl)
sample_ctrl[1,1:2,]

set.seed(905)
# random shuffle genes and samples
random_idx_gene = sample.int(nGeneTotal)
random_idx_sam  = sample.int(nall)
sample_ctrl = sample_ctrl[random_idx_gene, random_idx_sam, ]
RIN.simu = RIN.simu[random_idx_sam]

sample_param_case = sample_ctrl[, 1:ncase, ] 
sample_param_ctrl = sample_ctrl[, (ncase + 1):nall, ]

dim(sample_param_case)
dim(sample_param_ctrl)

sample_param_ctrl[1,1:2,]
sample_param_case[1,1:2,]

# check the mean value parameter across genes
rd_case = (colSums(sample_param_case[, , 1]))
rd_ctrl = (colSums(sample_param_ctrl[, , 1]))

t1 = t.test(rd_case, rd_ctrl)
t1

quantile(rd_case)
quantile(rd_ctrl)

# ------------------------------------------------------------------------
# modify the parameters to simulate differences between cases and controls
# ------------------------------------------------------------------------

# To make sure all parameters are non-negative, we do some transformation
r_mean2 = r_mean
r_var2  = r_var

if (r_mean > 1) {
  r_mean2 = 1 / r_mean
}

if (r_var < 1) {
  r_var2 = 1 / r_var
}

# ------------------------------------------------------------------------
# set up parameters for the genes with change in mean expression
# ------------------------------------------------------------------------

set.seed(2019)
runifs = rep(0, length(mean_index))
runifs[sample(length(mean_index), round(length(mean_index)/2))] = 1
runifs[1:9]
table(runifs)

# randomly choose half of the genes to modify in cases and 
# half of the genes to modify in controls
k = 0
for (i in mean_index) {
  k = k + 1
  if(runifs[k] > 0.5){
    for (j in 1:ncase) {
      sample_param_case[i, j, 1] = sample_param_case[i, j, 1]*r_mean2
    }
  }else{
    for (j in 1:nctrl) {
      sample_param_ctrl[i, j, 1] = sample_param_ctrl[i, j, 1]*r_mean2
    }
  }
}

mean_index[1:2]
runifs[1:2]
r_mean2
sample_param_case[mean_index[1], 1:5, ]
sample_param_ctrl[mean_index[1], 1:5, ]

rd_case    = colSums(sample_param_case[mean_index, , 1])
rd_control = colSums(sample_param_ctrl[mean_index, , 1])

summary(rd_case)
summary(rd_control)
t.test(rd_case, rd_control)

# ------------------------------------------------------------------------
# set up parameters for the genes with change in variance
# ------------------------------------------------------------------------

# the function calc_par_var returns the over-dispersion parameters 
# for changing the variance of a negative binomial distribution. 
# to make sure theta is larger than 0, r_v should be > 1. 

calc_par_var = function(mu, theta, r_v) {
  theta2 = theta  * mu / (mu * r_v + (r_v - 1) * theta)
  if(theta2 < 0){ stop("negative theta2.") }
  theta2
}

set.seed(2020)
runifs = rep(0, length(var_index))
runifs[sample(length(var_index), round(length(var_index)/2))] = 1
runifs[1:9]
table(runifs)

k = 0
for (i in var_index) {
  k = k + 1
  if(runifs[k] > 0.5){
    for (j in 1:ncase) {
      x = sample_param_case[i, j, ]
      sample_param_case[i,j,2] = calc_par_var(mu = x[1], theta = x[2],
                                              r_v = r_var2)
    }
  }else{
    for (j in 1:nctrl) {
      x = sample_param_ctrl[i, j, ]
      sample_param_ctrl[i,j,2] = calc_par_var(mu = x[1], theta = x[2],
                                              r_v = r_var2)
    }
  }
}

quantile(colSums(sample_param_case[var_index, , 1]))
quantile(colSums(sample_param_ctrl[var_index, , 1]))

# ------------------------------------------------------------------------
# check the parameters
# ------------------------------------------------------------------------

# check the total read depth
rd_case    = colSums(sample_param_case[, , 1])
rd_control = colSums(sample_param_ctrl[, , 1])

summary(rd_case)
summary(rd_control)
t.test(rd_case, rd_control)

# scatter plot
pdf(sprintf("figures/check_simulation_scatter_%s.pdf", config), 
    width = 9, height = 6)
par(mfrow = c(2, 3), pty = "s")

gene.list = list(non_DE=EE_index, mean_DE=mean_index, Var_DE=var_index)

for(k in 1:length(gene.list)){
  glistk = gene.list[[k]]
  nmk    = gsub("_", "-", names(gene.list)[k])
  plot(apply(log10(sample_param_ctrl[glistk,,1]), 1, median),
       apply(log10(sample_param_case[glistk,,1]), 1, median),
       cex = .2, xlab = "control cells", ylab = "case cells",
       main = sprintf("median of log10 mean, %s genes", nmk))
  abline(0, 1, col = "red")
}

for(k in 1:length(gene.list)){
  glistk = gene.list[[k]]
  nmk    = gsub("_", "-", names(gene.list)[k])
  plot(apply(log10(sample_param_ctrl[glistk,,2]), 1, median),
       apply(log10(sample_param_case[glistk,,2]), 1, median),
       cex = .2, xlab = "control cells", ylab = "case cells",
       main = sprintf("median of log10 dispersion, %s genes", nmk))
  abline(0, 1, col = "red")
}

dev.off()

# check the distribution of log_mean across genes, and sample_mean_sd_i

summary(log(sample_param_case[, 1, 1]))
summary(sample_param_case[, 1, 4])

summary(c(log(sample_param_case[, , 1])))
summary(c(sample_param_case[, , 4]))

# ------------------------------------------------------------------------
# simulate scRNAseq based on zinb parameters of cases and controls
# ------------------------------------------------------------------------
# Assume first ncase individuals are cases, the remaining are controls

if(UNEQ_N_CELL){
  ncell_each   = c(as.numeric(ncell_case), as.numeric(ncell_ctrl))
  ncell_all    = sum(ncell_each)
  ncell_cumsum = c(0, cumsum(ncell_each))
}else{
  ncell_each   = rep(ncell, nall)
  ncell_all    = sum(ncell_each)
  ncell_cumsum = c(0, cumsum(ncell_each))
}


sim_matrix = matrix(0, nrow = nGeneTotal, ncol = ncell_all)

set.seed(2018)
date()
for(i in 1:nall){
  if(i %% 5 ==0){
    cat(i, date(), "\n")
  }
  
  idx_i = (ncell_cumsum[i]+1):(ncell_cumsum[i+1])
  
  if (i > ncase) {
    mean_i = sample_param_ctrl[, (i - ncase), 1]
    disp_i = sample_param_ctrl[, (i - ncase), 2]
    drop_i = sample_param_ctrl[, (i - ncase), 3]
    sample_mean_sd_i = sample_param_ctrl[, (i - ncase), 4]
  } else{
    mean_i = sample_param_case[, i, 1]
    disp_i = sample_param_case[, i, 2]
    drop_i = sample_param_case[, i, 3]
    sample_mean_sd_i = sample_param_case[, i, 4]
  }
  
  sim_matrix[,idx_i] = 
  foreach(k = 1:ncell_each[i], .combine=cbind) %dorng% {
    sample_mean_k = exp(rnorm(nGeneTotal, log(mean_i), sample_mean_sd_i))
    sim_vector_cell_k = rep(NA, nGeneTotal)
    for (ig in 1:nGeneTotal) {
      sim_vector_cell_k[ig] = emdbook::rzinbinom(1, sample_mean_k[ig], 
                                                 disp_i[ig], drop_i[ig])
    }
    sim_vector_cell_k
  }
}
date()

dim(sim_matrix)
sim_matrix[1:8,1:6]

table(c(sim_matrix) == 0, useNA="ifany")
table(c(sim_matrix) == 0)/(nrow(sim_matrix)*ncol(sim_matrix))

####################### Meta information collection #################

# the phenotype and individual information of simulated samples.
phenotype  = c(rep(1, sum(ncell_each[1:ncase])), 
               rep(0, sum(ncell_each[(ncase+1):nall])))
individual = paste0("ind", c(rep(1:nall, times = ncell_each)))

table(phenotype)
sum(ncell_case)
sum(ncell_ctrl)

# Count info for matrix
cell_id = paste0("cell", 1:ncol(sim_matrix))
gene_id = paste0("gene", 1:nrow(sim_matrix))

rownames(sim_matrix) = gene_id
colnames(sim_matrix) = cell_id

# Cell info for meta
cell_rd = colSums(sim_matrix)
CDR     = colSums(sim_matrix > 0) / nrow(sim_matrix)
meta    = data.frame(cell_id, individual, phenotype, cell_rd, CDR, 
                     stringsAsFactors=FALSE)
dim(meta)
meta[1:2,]

meta_ind = meta[, c("individual", "phenotype")]
meta_ind = unique(meta_ind)
rownames(meta_ind) = meta_ind$individual

dim(meta_ind)
meta_ind[1:2,]

meta_ind$RIN = RIN.simu

pdf(sprintf("figures/check_covariates_%s.pdf", config), 
    width=6, height=3)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
boxplot(log10(meta$cell_rd) ~ meta$phenotype, xlab="group", 
        ylab="log10(read-depth)")
boxplot(meta$CDR ~ meta$phenotype, xlab="group", ylab="CDR")
dev.off()

# ------------------------------------------------------------------------
# save the simulated data
# ------------------------------------------------------------------------

dat_list = list(count_matrix = sim_matrix, meta_cell = meta, 
                meta_ind = meta_ind, gene_index = geneIndex)

saveRDS(dat_list, file=sprintf("data/sim_data_%s.rds", config))

sessionInfo()

mem_used()
gc()

q(save = "no")
