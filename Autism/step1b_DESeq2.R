
library(Matrix)
library(data.table)
library(dplyr)
library(DESeq2)


data.dir = "data"

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
dat1[1:5,1:4]

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
# generate individual level information
# ------------------------------------------------------------------------

meta_ind = distinct(meta[,3:12])
dim(meta_ind)
meta_ind[1:2,]
names(meta_ind)[9:10] = c("PMI", "RIN")

length(unique(meta$individual))

if(nrow(meta_ind) != length(unique(meta$individual))){
  stop("there is non-unique information\n")
}

table(meta_ind$Seqbatch, meta_ind$Capbatch)

# ------------------------------------------------------------------------
# collect count data
# ------------------------------------------------------------------------

trec1 = matrix(NA, nrow=nrow(dat1), ncol=nrow(meta_ind))
colnames(trec1) = meta_ind$sample
rownames(trec1) = rownames(dat1)
dim(trec1)
trec1[1:2,1:3]

for(i in 1:ncol(trec1)){
  wi = which(meta$sample == meta_ind$sample[i])
  trec1[,i] = rowSums(dat1[,wi])
}

dim(trec1)
trec1[1:2,1:3]

# ------------------------------------------------------------------------
# run DESeq2
# ------------------------------------------------------------------------

colData = meta_ind
for(i in 1:ncol(colData)){
  if(is.character(colData[[i]])){
    colData[[i]] = as.factor(colData[[i]])
  }
}
dim(colData)
colData[1:2,]
summary(colData)

colData$diagnosis = factor(colData$diagnosis, levels=c("Control", "ASD"))

dd0 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ diagnosis)
dd0  = DESeq(dd0)
res0 = results(dd0)
dim(res0)
head(res0)
summary(res0)


table(meta_ind$Seqbatch, meta_ind$Capbatch)
table(meta_ind$Seqbatch, meta_ind$diagnosis)

dds = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ age + sex + Seqbatch + RIN + diagnosis)
dds = DESeq(dds)

res = results(dds)
dim(res)
head(res)
summary(res)

nms = resultsNames(dds)
nms
nms = nms[-1]

pvals2 = matrix(NA, nrow=nrow(trec1), ncol=length(nms))

for(k in 1:length(nms)){
  rk = results(dds, name=nms[k])
  pvals2[,k] = rk$pvalue
}

colnames(pvals2) = nms
dim(pvals2)
head(pvals2)
summary(pvals2)

pdf(sprintf("figures/step1b_DESeq2_%s_pval_hist_final.pdf", grp), 
    width=9, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nms)){
  hist(pvals2[,k], main=nms[k], xlab="p-value", breaks=50)
}

dev.off()



# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(sprintf("figures/step1b_DESeq2_%s_compare_pval_Seq.pdf", grp), 
    width=9, height=3)
par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))
hist(res0$pvalue, main="without covariates", xlab="p-value")
hist(res$pvalue, main="with covariates (Seq)", xlab="p-value")
plot(-log10(res0$pvalue), -log10(res$pvalue), main="-log10(p-value)", 
     xlab="without covariates", ylab="with covariates (Seq)", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")
dev.off()



# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

dim(res0)
res0[1:2,]

dim(res)
res[1:2,]

res0 = as.data.frame(res0)
res  = as.data.frame(res)


fwrite(res0, file=sprintf("res/step1b_DESeq2_%s_no_covariates.tsv", grp), 
       sep = "\t")

fwrite(res, file=sprintf("res/step1b_DESeq2_%s_adj_covariates.tsv", grp), 
            sep = "\t")

gc()


sessionInfo()
q(save="no")
