# this code is modified from 
# 4f_neg_log_pvalue_vs_expression_level.R
# by using pvalue from PS_dca_direct_Was on a finer grid
# by setting number of permutations to 4999 and 99999

# the read pvalue part is also modified
# because we didn't rerun other parts for it, 
# so need to get pvalues directly from 1b DESeq2 and 5d dca_direct
# instead of getting them from the formal pvalue table 5e 
# we didn't construct this table because we don't have results
# from ideas and saver with the larger permutation number

# this code plots -log10(pvalue) against mean expression level
# any get a table for the values in PS_dca_direct_Was pvalues


# this file is updated from 
# step7h_update_plots_on_individual_gene_ideas_dca_direct
# under core_code/pipeline_formal folder

# the modifications:
# (1) deal with COVID data
# (2) add another way of adjusting for the read depth
#     following explanation in email


## current version handles CD8+Tcell 1 only


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

data.dir  = data.dir = "../../ideas_data/COVID/PBMC_10x"

grp = "CD8+Tcells_1"



# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"))
dim(cell_info)
cell_info[1:2,]

table(cell_info$id.celltype)

sort(table(paste(cell_info$group_per_sample, cell_info$donor, sep=":")))

# ------------------------------------------------------------------------
# read in count data of cell type
# ------------------------------------------------------------------------

dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat)
class(dat)
dat[1:5,1:4]

full_genes = rownames(dat) 

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat) %in% cell_info$cell)
meta = cell_info[match(colnames(dat), cell_info$cell),]
dim(meta)
meta[1:2,]


# ------------------------------------------------------------------------
# generate individual level information
# ------------------------------------------------------------------------


# filter out cells from control samples
meta_covid = meta[which(meta$group_per_sample != "control"),]
dim(meta_covid)

table(meta_covid$group_per_sample)
table(meta_covid$disease_stage)
table(meta_covid$donor)

table(meta_covid$donor, meta_covid$disease_stage)

df_donor = as.data.frame(table(meta_covid$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 10)]

meta2kp = meta_covid[which(meta_covid$donor %in% donor2kp),]
dim(meta2kp)
table(meta2kp$donor)
length(unique(meta2kp$donor))

cell2kp_index = which(meta$cell %in% meta2kp$cell)

# select counts of the cells to keep
dat1 = dat[, cell2kp_index]
mean(colnames(dat1) == meta2kp$cell)

meta_ind = distinct(meta2kp[,c('donor', 'group_per_sample', 'sex')])
table(meta_ind$group_per_sample, meta_ind$sex)


# ------------------------------------------------------------------------
# collect count data
# ------------------------------------------------------------------------

trec1 = matrix(NA, nrow=nrow(dat1), ncol=nrow(meta_ind))
colnames(trec1) = meta_ind$donor
rownames(trec1) = rownames(dat1)
dim(trec1)
trec1[1:2,1:3]

for(i in 1:ncol(trec1)){
  wi = which(meta2kp$donor == meta_ind$donor[i])
  trec1[,i] = rowSums(dat1[,wi])
}

dim(trec1)
trec1[1:2,1:3]



# find the row indexes to match those gene appearing in at least 90% of cells
n.zeros = rowSums(dat1 == 0)
summary(n.zeros)
w2kp = which(n.zeros < 0.9*ncol(dat1))
length(w2kp)




# DESeq2 pvalues did not filter out genes with large proportion of zeros
# need to repeat the filtering out genes step to get the genes we want
pvals_1b     = fread(sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_mild_severe.tsv", grp))
dim(pvals_1b)
pvals_1b$gene = full_genes
pvals_1b[1:2, ]

pvals_5d = fread(sprintf("res/5d_dca_direct_pvals_%s.tsv", grp))
dim(pvals_5d)
pvals_5d[1:2, ]

pvals_1b = pvals_1b[w2kp,]
dim(pvals_1b)

table(pvals_1b$gene == pvals_5d$gene)

table(full_genes[w2kp] == pvals_1b$gene)

pvals = data.frame(gene = pvals_1b$gene, 
                   DESeq2 = pvals_1b$pvalue, 
                   PS_dca_direct_Was = pvals_5d$PS_dca_direct_Was)


# check the distribution of small pvalues
df_dca_pvalue = as.data.frame(table(pvals$PS_dca_direct_Was))
colnames(df_dca_pvalue) = c("pvalue", "count")
df_dca_pvalue[1:6, ]

write.csv(df_dca_pvalue, 
          file = "res/5f_sorted_PS_dca_direct_Was_pvalue_table.csv", 
          row.names = FALSE)




# move on to get the rescaled counts based on the 
# description on page 14 of the DESeq2 paper
func_K_iR <- function(input_vec){
  if (min(input_vec) == 0){
    return(0)
  }
  else{
    log2_K_iR = mean(log2(input_vec))
    return(2^log2_K_iR)
  }
}

K_iR_vec = apply(trec1, 1, func_K_iR)
K_iR_vec = as.numeric(unlist(K_iR_vec))

s_j_vec = rep(1, ncol(trec1))
for (j in 1:ncol(trec1)){
  K_j_by_K_iR = as.numeric(unlist(trec1[, j]))/K_iR_vec
  # positive number divided by 0 will be Inf
  # 0 divided by 0 will be NA
  # remove the Inf and NAs to get the median
  K_j_by_K_iR_pure = K_j_by_K_iR[!is.na(K_j_by_K_iR) & !is.infinite(K_j_by_K_iR)]
  s_j_vec[j] = median(K_j_by_K_iR_pure)
}


scaled_trec1 = matrix(NA, nrow=nrow(trec1), ncol=nrow(meta_ind))
colnames(scaled_trec1) = colnames(trec1)
rownames(scaled_trec1) = rownames(trec1)
dim(scaled_trec1)
scaled_trec1[1:2, 1:3]

for (j in 1:ncol(trec1)){
  scaled_trec1[, j] = trec1[, j]/s_j_vec[j]
}



# subset to match the genes considered for ideas_nb and ideas_dca_direct
match_scaled_trec1 = scaled_trec1[w2kp, ]
match_scaled_trec1[1:2, ]

table(rownames(match_scaled_trec1) == pvals$gene)


mean_match_scaled_trec1 = apply(match_scaled_trec1, 1, mean)
mean_express_1 = as.numeric(mean_match_scaled_trec1)

# use the second way to do scaling to the trec1 matrix
djs = as.numeric(apply(trec1, 2, sum))
d_median = median(djs)
ajs = d_median/djs

scaled_trec2 = t(t(trec1) * ajs)
match_scaled_trec2 = scaled_trec2[w2kp, ]
table(rownames(match_scaled_trec2) == pvals$gene)

mean_match_scaled_trec2 = apply(match_scaled_trec2, 1, mean)
mean_express_2 = as.numeric(mean_match_scaled_trec2)


df_combined = data.frame(DESeq2 = -log10(pvals$DESeq2), 
                         dca_direct  = -log10(pvals$PS_dca_direct_Was), 
                         paper = log10(mean_express_1), 
                         email = log10(mean_express_2))

p1 = ggplot(df_combined, aes(x=email, y=DESeq2)) + 
  geom_point(alpha = 0.5, size = 1.2) + 
  xlab("log10(mean express) from email") + ylab("-log10(pval) DESeq2")

p2 = ggplot(df_combined, aes(x=email, y=dca_direct)) +
  geom_point(alpha = 0.5, size = 1.2) + 
  xlab("log10(mean express) from email") + ylab("-log10(pval) PS_dca_direct_Was") 

p3 = ggplot(df_combined, aes(x=paper, y=DESeq2)) + 
  geom_point(alpha = 0.5, size = 1.2) + 
  xlab("log10(mean express) from DESeq2") + ylab("-log10(pval) DESeq2")

p4 = ggplot(df_combined, aes(x=paper, y=dca_direct)) +
  geom_point(alpha = 0.5, size = 1.2) +  
  xlab("log10(mean express) from DESeq2") + ylab("-log10(pval) PS_dca_direct_Was") 


pdf("figures/5f_neglogpval_vs_mean_express.pdf", 
    width=6, height=6)
print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2))
dev.off()


df_scale_factors = data.frame("paper" = log10(s_j_vec),
                              "email" = -log10(ajs))
summary(df_scale_factors$paper - df_scale_factors$email)

p5 = ggplot(df_scale_factors, aes(x=paper, y=email)) +
  geom_point(alpha = 0.5, size = 1.2) +  
  xlab("log10(scale factor) from paper") + 
  ylab("log10(scale factor) from email") + 
  geom_abline()
p6 = ggplot(df_combined, aes(x=paper, y=email)) +
  geom_point(alpha = 0.5, size = 1.2) +  
  xlab("log10(mean express) from paper") + 
  ylab("log10(mean express) from email") + 
  geom_abline()

pdf("figures/5f_log_mean_express_email_vs_paper.pdf", 
    width=6, height=3)
print(ggarrange(p5, p6, ncol = 2, nrow = 1))
dev.off()


write.csv(df_combined, 
          file = "res/5f_mean_expressions.csv", 
          row.names = FALSE)


gc()

sessionInfo()
q(save="no")

