# compared to step1b_MAST.R, this version runs it on a subset to 
# get the computing time
# artificial way of doing subset

library(Matrix)
library(data.table)
library(dplyr)
library(MAST)
library(lme4)
library(doParallel)
library(foreach)
library(doRNG)

data.dir = "data"
nCore    = 4

args=(commandArgs(TRUE))
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
}else{
  eval(parse(text=args[[1]]))
}

grp

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"))
dim(cell_info)
cell_info[1:2,]

table(cell_info$region)

cell_info = cell_info[which(cell_info$region=="PFC"),]

sort(table(paste(cell_info$diagnosis, cell_info$sample, sep=":")))

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

dat1 = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat1)
class(dat1)

# get full gene
full_gene = row.names(dat1)

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat1) %in% cell_info$cell)
meta = cell_info[match(colnames(dat1), cell_info$cell),]
dim(meta)
meta[1:2,]

summary(meta$UMIs/meta$genes)

# check each individual has a unique sample
table(tapply(meta$sample, meta$individual, function(v){length(unique(v))}))

# check each individual has a unique Capbatch
table(tapply(meta$Capbatch, meta$individual, function(v){length(unique(v))}))

table(meta$cluster)
table(meta$region)

sort(table(paste(meta$diagnosis, meta$sample, sep=":")))

# ------------------------------------------------------------------------
# run MAST
# ------------------------------------------------------------------------

rds = colSums(dat1)
med_rds = median(rds)
summary(rds)
med_rds

dim(dat1)
dat1[1:3,1:6]
dat1 = t(t(dat1)/rds)*med_rds
dim(dat1)
dat1[1:3,1:6]
summary(colSums(dat1))

dat1_log = as.matrix(log2(1 + dat1)) #log transformed data

dim(dat1_log)
dat1_log[1:3, 1:4]
cell_id = colnames(dat1_log)   # get the cell id from the data
gene_id = rownames(dat1_log)   # get the gene id from the data



# do subset
df_load = fread("res/step1e_pvals_PFC_L2_3.tsv")
w2kp = match(df_load$gene, full_gene)


#fData = data.frame(primerid = gene_id)
fData = data.frame(primerid = gene_id[w2kp])
cData = data.frame(wellKey  = cell_id)

# Cell info for meta
cell_rd = colSums(dat1)
CDR     = colSums(dat1 > 0) / nrow(dat1)

meta_ind = meta[, c("individual", "diagnosis")]
meta_ind = unique(meta_ind)
rownames(meta_ind) = meta_ind$individual

dim(meta)
meta[1:2,]
dim(meta_ind)
meta_ind[1:2,]








gc()
#sca = FromMatrix(dat1_log, cData, fData)
sca = FromMatrix(dat1_log[w2kp,], cData, fData)
colData(sca)$cngeneson = as.numeric(CDR)
colData(sca)$diagnosis = as.factor(meta$diagnosis)
colData(sca)$ind = as.factor(meta$individual)
colData(sca)$RIN = meta$'RNA Integrity Number'
colData(sca)

#sca=sca[1:10,]

rm(dat1_log)
gc()

registerDoParallel(cores=nCore)
options(mc.cores=nCore)
getOption("mc.cores")

date()
b0 = zlm(formula = ~ diagnosis + cngeneson + RIN, sca = sca, 
         parallel = TRUE)
date()
b1 = zlm(formula = ~ diagnosis + (1 | ind) + cngeneson + RIN, sca = sca, 
         method = 'glmer', ebayes = FALSE, parallel = TRUE)
date()

b0
b1

date()
lrt0 = lrTest(b0, "diagnosis")
date()
lrt1 = lrTest(b1, "diagnosis")
date()

dim(lrt0)
lrt0[1,,]

dim(lrt1)
lrt1[1,,]

mast_pval_glm = apply(lrt0, 1, function(x){x[3,3]})
length(mast_pval_glm)
mast_pval_glm[1:4]

mast_pval_glmer = apply(lrt1, 1, function(x){x[3,3]})
length(mast_pval_glmer)
mast_pval_glmer[1:4]

# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(sprintf("figures/step1b_subset_MAST_%s_compare_pval_Seq.pdf", grp), 
    width=9, height=3)
par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))
hist(mast_pval_glm, main="MAST glm", xlab="p-value")
hist(mast_pval_glmer, main="MAST glmer", xlab="p-value")
plot(-log10(mast_pval_glm), -log10(mast_pval_glmer), main="-log10(p-value)", 
     xlab="MAST glm", ylab="MAST glmer", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")
dev.off()
# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------


fwrite(data.frame(mast_pval_glm), 
       file=sprintf("res/step1b_subset_MAST_%s_glm.tsv", grp), 
       sep = "\t", row.names=TRUE, col.names=FALSE)

fwrite(data.frame(mast_pval_glmer), 
       file=sprintf("res/step1b_subset_MAST_%s_glmer.tsv", grp), 
            sep = "\t", row.names=TRUE, col.names=FALSE)

gc()

sessionInfo()
q(save="no")
