
library(Matrix)
library(data.table)
library(dplyr)
library(DESeq2)


data.dir = "../../ideas_data/COVID/PBMC_10x"

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
# generate individual level information
# ------------------------------------------------------------------------

meta_ind_explore = distinct(meta[,c('donor', 'age', 'sex', 'group_per_sample')])
dim(meta_ind_explore)
# none of the control donors have either age or sex information
meta_ind_explore


length(unique(meta$donor))

if(nrow(meta_ind_explore) != length(unique(meta$donor))){
  stop("there is non-unique information\n")
}


# filter out cells from early stage samples

table(meta$donor, meta$disease_stage)

meta_no_early = meta[which(meta$disease_stage != "early"),]
dim(meta_no_early)

table(meta_no_early$donor, meta_no_early$disease_stage)

df_donor = as.data.frame(table(meta_no_early$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 30)]

meta2kp = meta_no_early[which(meta_no_early$donor %in% donor2kp),]
dim(meta2kp)

cell2kp_index = which(meta$cell %in% meta2kp$cell)

# select counts of cells to keep
dat1 = dat[, cell2kp_index]
mean(colnames(dat1) == meta2kp$cell)

meta_ind = distinct(meta2kp[,c('donor', 'group_per_sample')])
meta_ind$group_per_sample[which(meta_ind$group_per_sample != "control")] = "COVID"
meta_ind

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


# ------------------------------------------------------------------------
# run DESeq2
# ------------------------------------------------------------------------

colData = meta_ind
colnames(colData)[2] = 'diagnosis'
for(i in 1:ncol(colData)){
  if(is.character(colData[[i]])){
    colData[[i]] = as.factor(colData[[i]])
  }
}
dim(colData)
colData[1:2,]


summary(colData)

colData$diagnosis = factor(colData$diagnosis, levels=c("control", "COVID"))


# ------------------------------------------------------------------------
# add to colData individual level covariates to adjust for 
# ------------------------------------------------------------------------

# total cell number under each individual

n_cells = c()

for (i in c(1:dim(meta_ind)[1])){
  n_cells = c(n_cells, length(which(meta2kp$donor == meta_ind$donor[i])))
}

colData$ncells = n_cells

# average cell level read depth under each individual

average_cell_rd = c()

for (i in c(1:dim(meta_ind)[1])){
  wi = which(meta2kp$donor == meta_ind$donor[i])
  cell_rds = apply(dat1[,wi], 2, sum)
  average_cell_rd = c(average_cell_rd, mean(cell_rds))
}

colData$rd = average_cell_rd







dd0 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ diagnosis)
dd0  = DESeq(dd0)
res0 = results(dd0)
dim(res0)
head(res0)
summary(res0)

nm0 = resultsNames(dd0)
nm0
nm0 = nm0[-1]

pvals0 = matrix(NA, nrow=nrow(trec1), ncol=length(nm0))

for(k in 1:length(nm0)){
  rk = results(dd0, name=nm0[k])
  pvals0[,k] = rk$pvalue
}

colnames(pvals0) = nm0
dim(pvals0)
head(pvals0)
summary(pvals0)

pdf(sprintf("figures/1b_DESeq2_%s_pval_hist_no_covariate_repeat.pdf", grp), 
    width=4.5, height=3.5)
hist(pvals0[,1], main=nm0[1], xlab="p-value", breaks=50)
dev.off()




# first, only include average cell level read depth as covariate

dd1 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ log(rd) + diagnosis)
dd1 = DESeq(dd1)

res1 = results(dd1)
dim(res1)
head(res1)
summary(res1)

nm1 = resultsNames(dd1)
nm1
nm1 = nm1[-1]

pvals1 = matrix(NA, nrow=nrow(trec1), ncol=length(nm1))

for(k in 1:length(nm1)){
  rk = results(dd1, name=nm1[k])
  pvals1[,k] = rk$pvalue
}

colnames(pvals1) = nm1
dim(pvals1)
head(pvals1)
summary(pvals1)


pdf(sprintf("figures/1b_DESeq2_%s_pval_hist_logrd.pdf", grp), 
    width=9, height=3.5)
par(mfrow=c(1,2), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nm1)){
  hist(pvals1[,k], main=nm1[k], xlab="p-value", breaks=50)
}

dev.off()



# second, include both number of cells and 
# average cell level read depth as covariates

dd2 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ log(ncells) + log(rd) + diagnosis)
dd2 = DESeq(dd2)

res2 = results(dd2)
dim(res2)
head(res2)
summary(res2)

nm2 = resultsNames(dd2)
nm2
nm2 = nm2[-1]

pvals2 = matrix(NA, nrow=nrow(trec1), ncol=length(nm2))

for(k in 1:length(nm2)){
  rk = results(dd2, name=nm2[k])
  pvals2[,k] = rk$pvalue
}

colnames(pvals2) = nm2
dim(pvals2)
head(pvals2)
summary(pvals2)


pdf(sprintf("figures/1b_DESeq2_%s_pval_hist_logncells_logrd.pdf", grp), 
    width=13.5, height=3.5)
par(mfrow=c(1,3), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nm2)){
  hist(pvals2[,k], main=nm2[k], xlab="p-value", breaks=50)
}

dev.off()



# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(sprintf("figures/1b_DESeq2_%s_compare_pval.pdf", grp), 
    width=9, height=3)
par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))

plot(-log10(res0$pvalue), -log10(res1$pvalue), main="-log10(p-value)", 
     xlab="without covariates", ylab="with log(rd)", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

plot(-log10(res0$pvalue), -log10(res2$pvalue), main="-log10(p-value)", 
     xlab="without covariates", ylab="with log(ncells)+log(rd)", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

plot(-log10(res1$pvalue), -log10(res2$pvalue), main="-log10(p-value)", 
     xlab="with log(rd)", ylab="with log(ncells)+log(rd)", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

dev.off()



# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

dim(res0)
res0[1:2,]

dim(res1)
res1[1:2,]

dim(res2)
res2[1:2,]

res0 = as.data.frame(res0)
res1 = as.data.frame(res1)
res2 = as.data.frame(res2)

fwrite(res0, file=sprintf("res/1b_DESeq2_%s_no_covariates_repeat.tsv", grp), 
       sep = "\t")

fwrite(res1, file=sprintf("res/1b_DESeq2_%s_logrd.tsv", grp), 
            sep = "\t")

fwrite(res2, file=sprintf("res/1b_DESeq2_%s_logncells_logrd.tsv", grp), 
       sep = "\t")

gc()


sessionInfo()
q(save="no")
