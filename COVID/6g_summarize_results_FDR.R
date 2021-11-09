# this file gets the number of significant genes identified by
# different methods at certain FDR cutoffs


library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(stringr)
library(tidyr)


data.dir = "../../ideas_data/COVID/PBMC_10x"


# cell type information
cell_types = c("CD8+Tcells_1")
grp = cell_types[1]

# method list
methods = c("rank-sum", "MAST", "MAST-glmer", "DESeq2", "IDEAS", "DCA", "SAVER")

# get the genes to keep
# read in cell information
cell_info = fread(file.path(data.dir, "meta.tsv"), 
                  stringsAsFactors=TRUE)
dim(cell_info)
cell_info[1:2,]

# read in count data of one celltype
dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat)
class(dat)
dat[1:5,1:4]

# read in covid donor information
covid_donor_info = 
  read.csv(file.path(data.dir, "covid_donor_info_from_mmc1.csv"), 
           header = TRUE)

dim(covid_donor_info)
covid_donor_info[1:2,]
summary(covid_donor_info)

# subset cell information
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

# filter out genes with too many zero's
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


# -------------------------------------------------------------------
# load results of different methods after permutations
# -------------------------------------------------------------------



pv_perm = list()

DESeq2     = fread(sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_perm.tsv", grp))
mast_glm   = fread(sprintf("res/3b_MAST_perm_%s_glm.tsv", grp))
mast_glmer = fread(sprintf("res/3b_MAST_perm_%s_glmer.tsv", grp))
rank_sum   = fread(sprintf("res/1b_ranksum_perm_%s.tsv", grp))
saver      = fread(sprintf("res/5d_saver_direct_pvals_%s_perm.tsv", grp))
ideas      = fread(sprintf("res/5c_nb_pvals_%s_perm.tsv", grp))
dca        = fread(sprintf("res/5d_dca_direct_pvals_%s_perm.tsv", grp))

DESeq2 = DESeq2[w2kp,]
rank_sum   = rank_sum[w2kp,]

table(rownames(dat1) == mast_glm$V1)
table(rownames(dat1) == mast_glmer$V1) 
table(rownames(dat1) == ideas$gene) 
table(rownames(dat1) == dca$gene) 
table(rownames(dat1) == saver$gene)

pv_perm[["rank-sum"]]   = rank_sum$V2
pv_perm[["MAST"]]       = mast_glm$V2
pv_perm[["MAST-glmer"]] = mast_glmer$V2
pv_perm[["DESeq2"]]     = DESeq2$pvalue
pv_perm[["IDEAS"]]      = ideas$PS_nb_Was
pv_perm[["DCA"]]        = dca$PS_dca_direct_Was
pv_perm[["SAVER"]]      = saver$PS_saver_direct_Was




sorted_pv_perm = list()

for (m in methods){
  sorted_pv_perm[[m]] = sort(pv_perm[[m]])
}


get_num_NA <- function(vec){
  return(sum(is.na(vec)))
}

sapply(pv_perm, get_num_NA)
sapply(sorted_pv_perm, length)

# -------------------------------------------------------------------
# load true pvalue results from different methods
# -------------------------------------------------------------------

p1  = fread(sprintf("res/5e_pvals_%s.tsv", grp))

pv = c()

pv[["rank-sum"]]   = p1$rank_sum
pv[["MAST"]]       = p1$MAST_glm
pv[["MAST-glmer"]] = p1$MAST_glmer
pv[["DESeq2"]]     = p1$DESeq2
pv[["IDEAS"]]      = p1$PS_nb_Was
pv[["DCA"]]        = p1$PS_dca_direct_Was
pv[["SAVER"]]      = p1$PS_saver_direct_Was


sapply(pv, get_num_NA)
sapply(pv, length)


sorted_pv = list()

for (m in methods){
  sorted_pv[[m]] = sort(pv[[m]])
}


sapply(sorted_pv, get_num_NA)
sapply(sorted_pv, length)

# -------------------------------------------------------------------
# count number of false discoveries at each cutoffs
# -------------------------------------------------------------------

cnt = list()

for (m in methods){
  cur_sorted_pvalues = sorted_pv[[m]]
  cnt_fds = rep(NA, length(cur_sorted_pvalues))
  for (i in 1:length(cur_sorted_pvalues)){
    cnt_fds[i] = sum(sorted_pv_perm[[m]] <= cur_sorted_pvalues[i])
  }
  cnt[[m]] = cnt_fds
}


# -------------------------------------------------------------------
# compute q hat
# -------------------------------------------------------------------

q_hat = list()

for (m in methods){
  n = length(sorted_pv[[m]])
  q_hat_vec = rep(NA, n)
  q_hat_vec[n] = cnt[[m]][n]/n
  for (j in 1:(n-1)){
    k = n - j
    q_hat_vec[k] = min(cnt[[m]][k]/k, q_hat_vec[k+1])
  }
  q_hat[[m]] = q_hat_vec
}


saveRDS(q_hat, file = "res/6g_q_hat_list.rds")

sapply(q_hat, function(x){sum(x < 0.05)})


# -------------------------------------------------------------------
# get number of significant genes from different methods
# under different FDR cutoffs
# -------------------------------------------------------------------

FDR_cuts = c(0.001, 0.005, 0.01, 0.05, 0.1)

sign_n_mat = matrix(NA, 
                    nrow = length(methods), 
                    ncol = length(FDR_cuts))
for (j in 1:length(FDR_cuts)){
  sign_n_mat[, j] = sapply(q_hat, function(x){sum(x < FDR_cuts[j])})
}

df_sign_num = as.data.frame(cbind(methods, sign_n_mat))
colnames(df_sign_num)[2:ncol(df_sign_num)] = FDR_cuts

write.csv(df_sign_num, 
          file = "res/6g_num_sign_under_FDR_cutffs.csv", 
          row.names = FALSE)



gc()

sessionInfo()
q(save="no")
