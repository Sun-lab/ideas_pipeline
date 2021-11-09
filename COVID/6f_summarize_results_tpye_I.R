
library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(stringr)
library(tidyr)

theme_set(theme_classic())

data.dir = "../../ideas_data/COVID/PBMC_10x"


# -------------------------------------------------------------------
# cell type information
# -------------------------------------------------------------------
cell_types = c("CD8+Tcells_1")

# -------------------------------------------------------------------
# load MAST, rank-sum, and ideas results after permutations
# -------------------------------------------------------------------

pcuts = c(1e-5, 1e-4, 0.001, 0.01, 0.05)

type_i_mast_glm = matrix(NA, nrow=length(cell_types), ncol=length(pcuts))
rownames(type_i_mast_glm) = cell_types
colnames(type_i_mast_glm) = pcuts
type_i_mast_glm

type_i = list()

methods = c("rank-sum", "MAST", "MAST-glmer", "DESeq2", "IDEAS", "DCA", "SAVER")
for(m1 in methods){
  type_i[[m1]] = type_i_mast_glm
}


for(k in 1:length(cell_types)){
  
  grp = cell_types[k]
  
  # ---------------
  # first part is to get the genes to keep
  # ---------------
  
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
  
  # ----------------
  # the second part loads the pvalues from permutation
  # ----------------
  
  # load permutation pvalues
  
  DESeq2     = fread(sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_perm.tsv", grp))
  mast_glm   = fread(sprintf("res/3b_MAST_perm_%s_glm.tsv", grp))
  mast_glmer = fread(sprintf("res/3b_MAST_perm_%s_glmer.tsv", grp))
  rank_sum   = fread(sprintf("res/1b_ranksum_perm_%s.tsv", grp))
  saver      = fread(sprintf("res/5d_saver_direct_pvals_%s_perm.tsv", grp))
  ideas      = fread(sprintf("res/5c_nb_pvals_%s_perm.tsv", grp))
  dca        = fread(sprintf("res/5d_dca_direct_pvals_%s_perm.tsv", grp))
  
  dim(ideas)
  ideas[1:2,]
  
  table(rownames(dat1) == mast_glm$V1)
  table(rownames(dat1) == mast_glmer$V1) 
  table(rownames(dat1) == ideas$gene) 
  table(rownames(dat1) == dca$gene) 
  table(rownames(dat1) == saver$gene)

  
  DESeq2   = DESeq2[w2kp,]
  rank_sum = rank_sum[w2kp,]
  
  for(j in 1:length(pcuts)){
    pj = pcuts[j]
    type_i$`rank-sum`[k,j]   = mean(rank_sum$V2 < pcuts[j], na.rm = TRUE)
    type_i$MAST[k,j]         = mean(mast_glm$V2 < pcuts[j], na.rm = TRUE)
    type_i$`MAST-glmer`[k,j] = mean(mast_glmer$V2 < pcuts[j], na.rm = TRUE)
    type_i$DESeq2[k,j]       = mean(DESeq2$pvalue < pcuts[j], na.rm = TRUE)
    type_i$IDEAS[k,j]        = mean(ideas$PS_nb_Was < pcuts[j], na.rm = TRUE)
    type_i$DCA[k,j]          = mean(dca$PS_dca_direct_Was < pcuts[j], na.rm = TRUE)
    type_i$SAVER[k,j]        = mean(saver$PS_saver_direct_Was < pcuts[j], na.rm = TRUE)
  }
}

lapply(type_i, dim)
type_i


# save type_I error results out

type_i_mat = matrix(NA, nrow = 7, ncol = 5)
for (i in 1:length(methods)){
  m = methods[i]
  type_i_mat[i, ] = type_i[[m]][1, ]
}

df_type_i = 
  as.data.frame(cbind(methods, signif(type_i_mat, digits = 2)))
colnames(df_type_i)[2:ncol(df_type_i)] = pcuts

write.csv(df_type_i, file = "res/6f_type_i.csv", 
          row.names = FALSE)

gc()

sessionInfo()
q(save="no")
