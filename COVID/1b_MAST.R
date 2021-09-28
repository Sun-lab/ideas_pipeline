
library(Matrix)
library(data.table)
library(dplyr)
library(MAST)
library(lme4)
library(doParallel)
library(foreach)
library(doRNG)

data.dir = "../../ideas_data/COVID/PBMC_10x"
nCore    = 4

args=(commandArgs(TRUE))
args

if (length(args) != 1) {
  message("one argument is expected, use 'CD8+Tcells_1' as default.\n")
  grp = "CD8+Tcells_1"
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

table(cell_info$id.celltype)

sort(table(paste(cell_info$group_per_sample, cell_info$donor, sep=":")))

# ------------------------------------------------------------------------
# read in count data of cell type
# ------------------------------------------------------------------------

dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat)
class(dat)
dat[1:5,1:4]

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat) %in% cell_info$cell)
meta = cell_info[match(colnames(dat), cell_info$cell),]
dim(meta)
meta[1:2,]

summary(meta$nCount_RNA/meta$nFeature_RNA)

# check how many samples (each sample contains multiple cells) each donor has
table(tapply(meta$sampleID, meta$donor, function(v){length(unique(v))}))


table(meta$cluster_labels_res.0.4)
table(meta$id.celltype)

sort(table(paste(meta$group_per_sample, meta$sampleID, sep=":")))


# ------------------------------------------------------------------------
# filter out cells from early stage samples
# get individual level information
# ------------------------------------------------------------------------

table(meta$donor, meta$disease_stage)

meta_no_early = meta[which(meta$disease_stage != "early"),]
dim(meta_no_early)

table(meta_no_early$donor, meta_no_early$disease_stage)

df_donor = as.data.frame(table(meta_no_early$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 30)]

meta2kp = meta_no_early[which(meta_no_early$donor %in% donor2kp),]
dim(meta2kp)

# create a column for cell level label in terms of COVID and control
table(meta2kp$group_per_sample)
meta2kp$diagnosis = meta2kp$group_per_sample
meta2kp$diagnosis[which(meta2kp$group_per_sample == "mild")] = "COVID"
meta2kp$diagnosis[which(meta2kp$group_per_sample == "severe")] = "COVID"
table(meta2kp$diagnosis)
dim(meta2kp)

cell2kp_index = which(meta$cell %in% meta2kp$cell)

# select counts of cells to keep
dat1 = dat[, cell2kp_index]
mean(colnames(dat1) == meta2kp$cell)

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

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey  = cell_id)

# Cell info for meta
cell_rd = colSums(dat1)
CDR     = colSums(dat1 > 0) / nrow(dat1)


meta_ind = distinct(meta2kp[,c('donor', 'group_per_sample')])
meta_ind$group_per_sample[which(meta_ind$group_per_sample != "control")] = "COVID"
colnames(meta_ind)[2] = 'diagnosis'
rownames(meta_ind) = meta_ind$donor
meta_ind


dim(meta2kp)
meta2kp[1:2,]
dim(meta_ind)
meta_ind[1:2,]

gc()
sca = FromMatrix(dat1_log, cData, fData)
colData(sca)$cngeneson = as.numeric(CDR)
colData(sca)$diagnosis = as.factor(meta2kp$diagnosis)
colData(sca)$ind = as.factor(meta2kp$donor)
#colData(sca)$RIN = meta$'RNA Integrity Number'
colData(sca)

#sca=sca[1:10,]

rm(dat1_log)
gc()

registerDoParallel(cores=nCore)
options(mc.cores=nCore)
getOption("mc.cores")

date()
b0 = zlm(formula = ~ diagnosis + cngeneson, sca = sca, 
         parallel = TRUE)
#b0 = zlm(formula = ~ diagnosis + cngeneson + RIN, sca = sca, 
#         parallel = TRUE)
date()
b1 = zlm(formula = ~ diagnosis + (1 | ind) + cngeneson, sca = sca, 
         method = 'glmer', ebayes = FALSE, parallel = TRUE)
#b1 = zlm(formula = ~ diagnosis + (1 | ind) + cngeneson + RIN, sca = sca, 
#         method = 'glmer', ebayes = FALSE, parallel = TRUE)
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

pdf(sprintf("figures/1b_MAST_%s_compare_pval_Seq.pdf", grp), 
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
       file=sprintf("res/1b_MAST_%s_glm.tsv", grp), 
       sep = "\t", row.names=TRUE, col.names=FALSE)

fwrite(data.frame(mast_pval_glmer), 
       file=sprintf("res/1b_MAST_%s_glmer.tsv", grp), 
            sep = "\t", row.names=TRUE, col.names=FALSE)

gc()

sessionInfo()
q(save="no")
