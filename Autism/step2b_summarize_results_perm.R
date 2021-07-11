
library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())
data.dir  = "./data"



qqp <- function(pvals, main, confidence=.95, cutoff=1){
  
  alpha = 1-confidence
  n     = length(pvals)
  
  pvals[is.na(pvals)]=1
  pvals=sort(pvals)
  
  k=c(1:n)
  
  lower = cutoff*qbeta(alpha/2, k, n+1-k)
  upper = cutoff*qbeta((1-alpha/2), k, n+1-k)
  
  expected = cutoff*k/(n+1)
  n0 = length(which(pvals ==0))
  
  if(n0 > 0){
    warning(sprintf("there are %d p-values being 0\n", n0))
  }
  
  biggest= max(-log10(pvals[which(pvals > 0)]), -log10(expected))
  
  plot(-log10(expected), -log10(pvals), xlim=c(0,biggest),
       ylim=c(0,biggest), pch=20, xlab="-log10(expected p-value)",
       ylab="-log10(observed p-value)", cex=0.6, bty="n", main=main)
  
  lines(-log10(expected), -log10(lower), lty=2)
  lines(-log10(expected), -log10(upper), lty=2)
  
}

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

# -------------------------------------------------------------------
# load MAST and rank sum results after permutations
# -------------------------------------------------------------------

pcuts = c(1e-5, 1e-4, 0.001, 0.01, 0.05)

type_i_mast_glm = matrix(NA, nrow=length(cell_types), ncol=length(pcuts))
rownames(type_i_mast_glm) = cell_types
colnames(type_i_mast_glm) = pcuts
type_i_mast_glm[1:2,]

type_i_mast_glmer = type_i_mast_glm
type_i_rank_sum   = type_i_mast_glm

for(k in 1:length(cell_types)){
  grp = cell_types[k]
  mast_glm   = fread(sprintf("res/step1b_MAST_perm_%s_glm.tsv", grp))
  mast_glmer = fread(sprintf("res/step1b_MAST_perm_%s_glmer.tsv", grp))
  rank_sum   = fread(sprintf("res/step1b_ranksum_perm_%s.tsv", grp))

  full_dat1  = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
  full_genes = rownames(full_dat1)  
  
  n.zeros = rowSums(full_dat1 == 0)
  summary(n.zeros)
  w2kp = which(n.zeros < 0.8*ncol(full_dat1))
  length(w2kp)
  
  stopifnot(all(full_genes == mast_glm$V1))
  stopifnot(all(full_genes == mast_glmer$V1))
  
  mast_glm   = mast_glm[w2kp,]
  mast_glmer = mast_glmer[w2kp,]
  rank_sum   = rank_sum[w2kp,]
  rank_sum_v = unlist(rank_sum[,-1])
  
  for(j in 1:length(pcuts)){
    pj = pcuts[j]
    type_i_mast_glm[k,j]   = mean(mast_glm$V2 < pcuts[j], na.rm = TRUE)
    type_i_mast_glmer[k,j] = mean(mast_glmer$V2 < pcuts[j], na.rm = TRUE)
    type_i_rank_sum[k,j]   = mean(rank_sum_v < pcuts[j], na.rm = TRUE)
  }
  
}

dim(type_i_mast_glm)
dim(type_i_mast_glmer)
dim(type_i_rank_sum)

type_i_mast_glm
type_i_mast_glmer
type_i_rank_sum

# -------------------------------------------------------------------
# load rank-sum results after permutations
# -------------------------------------------------------------------

gc()

sessionInfo()
q(save="no")
