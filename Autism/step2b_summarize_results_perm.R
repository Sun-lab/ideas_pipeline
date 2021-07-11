
library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(stringr)

theme_set(theme_classic())
data.dir  = "./data"

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

# -------------------------------------------------------------------
# load MAST, rank-sum, and ideas results after permutations
# -------------------------------------------------------------------

pcuts = c(1e-5, 1e-4, 0.001, 0.01, 0.05)

type_i_mast_glm = matrix(NA, nrow=length(cell_types), ncol=length(pcuts))
rownames(type_i_mast_glm) = cell_types
colnames(type_i_mast_glm) = pcuts
type_i_mast_glm[1:2,]

type_i = list()

methods = c("rank-sum", "MAST", "MAST-glmer", "DESeq2", "IDEAS", "DCA", "SAVER")
for(m1 in methods){
  type_i[[m1]] = type_i_mast_glm
}

for(k in 1:length(cell_types)){
  grp = cell_types[k]
  mast_glm   = fread(sprintf("res/step1b_MAST_perm_%s_glm.tsv", grp))
  mast_glmer = fread(sprintf("res/step1b_MAST_perm_%s_glmer.tsv", grp))
  rank_sum   = fread(sprintf("res/step1b_ranksum_perm_%s.tsv", grp))
  saver      = fread(sprintf("res/step1d_saver_direct_perm_pvals_%s.tsv", grp))
  ideas      = fread(sprintf("permutations/res/step11r_pvals_%s_ten_p.tsv", grp))
  
  dim(ideas)
  ideas[1:2,]
  
  full_dat1  = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
  dn         = dimnames(full_dat1)
  full_genes = dn[[1]]
  
  n.zeros = rowSums(full_dat1 == 0)
  summary(n.zeros)
  w2kp = which(n.zeros < 0.8*ncol(full_dat1))
  length(w2kp)
  
  stopifnot(all(full_genes == mast_glm$V1))
  stopifnot(all(full_genes == mast_glmer$V1))
  stopifnot(all(full_genes[w2kp] == saver$gene))
  # somewhat the dimnames or rownames does not work within a loop
  # stopifnot(setequal(full_genes[w2kp], ideas$gene))
  
  mast_glm   = mast_glm[w2kp,]
  mast_glmer = mast_glmer[w2kp,]
  rank_sum   = rank_sum[w2kp,]
  rank_sum_v = unlist(rank_sum[,-1])
  
  for(j in 1:length(pcuts)){
    pj = pcuts[j]
    type_i$`rank-sum`[k,j]   = mean(rank_sum_v < pcuts[j], na.rm = TRUE)
    type_i$MAST[k,j]         = mean(mast_glm$V2 < pcuts[j], na.rm = TRUE)
    type_i$`MAST-glmer`[k,j] = mean(mast_glmer$V2 < pcuts[j], na.rm = TRUE)
    type_i$DESeq2[k,j]       = mean(ideas$DESeq2 < pcuts[j], na.rm = TRUE)
    type_i$IDEAS[k,j]        = mean(ideas$PS_nb_Was < pcuts[j], na.rm = TRUE)
    type_i$DCA[k,j]          = mean(ideas$PS_dca_direct_Was < pcuts[j], na.rm = TRUE)
    type_i$SAVER[k,j]        = mean(saver$PS_saver_direct_Was < pcuts[j], na.rm = TRUE)
  }
}

lapply(type_i, dim)
type_i

# -------------------------------------------------------------------
# plot type I error across cell types
# -------------------------------------------------------------------

idx = 5
df1 = lapply(type_i, function(x){x[,idx]})
df1 = unlist(df1)
df2 = data.frame(method = str_extract(names(df1), "\\S+(?=.PFC_)"),
                 cell_type = str_extract(names(df1), "(?<=.PFC_)\\S+"),
                 type_I_error = as.numeric(df1))

df2$method = factor(df2$method, levels=methods)

p5 = ggplot(df2, aes(x=method, y=type_I_error, fill=method)) +
  geom_boxplot() + geom_hline(yintercept=0.05) + 
  geom_hline(yintercept=0.01, linetype="dashed") + 
  geom_hline(yintercept=0.10, linetype="dashed") + 
  scale_y_continuous(trans='log10') + ylab("Type-I error") + 
  theme(axis.text.x = element_blank())

idx = 4
df1 = lapply(type_i, function(x){x[,idx]})
df1 = unlist(df1)
df2 = data.frame(method = str_extract(names(df1), "\\S+(?=.PFC_)"),
                 cell_type = str_extract(names(df1), "(?<=.PFC_)\\S+"),
                 type_I_error = as.numeric(df1))

df2$method = factor(df2$method, levels=methods)
summary(df2$type_I_error)
min(df2$type_I_error[df2$type_I_error > 0])
df2$type_I_error[which(df2$type_I_error == 0)] = 0.0002

p1 = ggplot(df2, aes(x=method, y=type_I_error, fill=method)) +
  geom_boxplot() + geom_hline(yintercept=0.01) + 
  geom_hline(yintercept=0.001, linetype="dashed") + 
  geom_hline(yintercept=0.05, linetype="dashed") + 
  scale_y_continuous(trans='log10') + ylab("Type-I error") + 
  theme(axis.text.x  = element_blank())

pdf("figures/step2b_type_i_error_05.pdf", width=4, height=3)
p5
dev.off()

pdf("figures/step2b_type_i_error_01.pdf", width=4, height=3)
p1
dev.off()


gc()

sessionInfo()
q(save="no")
