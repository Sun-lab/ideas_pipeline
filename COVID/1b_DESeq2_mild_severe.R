
library(Matrix)
library(data.table)
library(dplyr)
library(DESeq2)


data.dir = "../ideas_data/COVID/PBMC_10x"

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
# read in covid donor information
# ------------------------------------------------------------------------

covid_donor_info = 
  read.csv(file.path(data.dir, "covid_donor_info_from_mmc1.csv"), 
           header = TRUE)

dim(covid_donor_info)
covid_donor_info[1:2,]
summary(covid_donor_info)

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



# check whether the control label from disease_stage and that from
# group_per_sample columns match
table(meta$disease_stage)
table(meta$group_per_sample)
table(meta$disease_stage, meta$group_per_sample)

table(meta$donor, meta$disease_stage)

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

# add exact age information
donor_info_match = 
  covid_donor_info[match(meta_ind$donor, covid_donor_info$donor),]
# double check that the condition and sex features match
mean(donor_info_match$condition == meta_ind$group_per_sample)
mean(which(donor_info_match$sex=="f") == which(meta_ind$sex=="female"))
meta_ind$age = donor_info_match$age

sort(meta_ind$age[which(meta_ind$group_per_sample=="mild")])
sort(meta_ind$age[which(meta_ind$group_per_sample=="severe")])

table(meta_ind$group_per_sample)

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
# keep the age column being numeric
for(i in 1:(ncol(colData)-1)){
  if(is.character(colData[[i]])){
    colData[[i]] = as.factor(colData[[i]])
  }
}
dim(colData)
colData[1:2,]


summary(colData)

colData$diagnosis = factor(colData$diagnosis, levels=c("mild", "severe"))


# ------------------------------------------------------------------------
# add to colData donor level covariates to adjust for 
# ------------------------------------------------------------------------

# total read depth across all cells under each individual

total_rd = c()

for (i in c(1:dim(meta_ind)[1])){
  wi = which(meta2kp$donor == meta_ind$donor[i])
  total_rd = c(total_rd, sum(apply(dat1[,wi], 2, sum)))
}

colData$totalrd = total_rd



# first, null model

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



# second, only include donor level total read depth as covariate

dd1 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ log(totalrd) + diagnosis)
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


# third, include sex and age as covariates

dd2 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ sex + age + diagnosis)
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



# fourth, include both donor level total read depth, sex and age

dd3 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ log(totalrd) + sex + age + diagnosis)
dd3 = DESeq(dd3)

res3 = results(dd3)
dim(res3)
head(res3)
summary(res3)

nm3 = resultsNames(dd3)
nm3
nm3 = nm3[-1]

pvals3 = matrix(NA, nrow=nrow(trec1), ncol=length(nm3))

for(k in 1:length(nm3)){
  rk = results(dd3, name=nm3[k])
  pvals3[,k] = rk$pvalue
}

colnames(pvals3) = nm3
dim(pvals3)
head(pvals3)
summary(pvals3)




pdf(sprintf("figures/1b_DESeq2_%s_pval_hists_mild_severe.pdf", grp), 
    width=18, height=3.5)
par(mfrow=c(1,4), bty="n", mar=c(5,4,2,1))
k = length(nm0)
hist(pvals0[,k], main="null", xlab="p-value", breaks=50)
k = length(nm1)
hist(pvals1[,k], main="log(totalrd)", xlab="p-value", breaks=50)
k = length(nm2)
hist(pvals2[,k], main="sex+age", xlab="p-value", breaks=50)
k = length(nm3)
hist(pvals3[,k], main="log(totalrd)+sex+age", xlab="p-value", breaks=50)
dev.off()



# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(sprintf("figures/1b_DESeq2_%s_compare_pval_mild_severe.pdf", grp), 
    width=9, height=3)
par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))

plot(-log10(res0$pvalue), -log10(res1$pvalue), main="-log10(p-value)", 
     xlab="null", ylab="log(totalrd)", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

plot(-log10(res0$pvalue), -log10(res2$pvalue), main="-log10(p-value)", 
     xlab="null", ylab="sex+age", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

plot(-log10(res0$pvalue), -log10(res3$pvalue), main="-log10(p-value)", 
     xlab="null", ylab="log(totalrd)+sex+age", 
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

dim(res3)
res3[1:2,]

res0 = as.data.frame(res0)
res1 = as.data.frame(res1)
res2 = as.data.frame(res2)
res3 = as.data.frame(res3)


fwrite(res0, file=sprintf("res/1b_DESeq2_%s_no_covariates_mild_severe.tsv", grp), 
       sep = "\t")

fwrite(res1, file=sprintf("res/1b_DESeq2_%s_logtotalrd_mild_severe.tsv", grp), 
       sep = "\t")

fwrite(res2, file=sprintf("res/1b_DESeq2_%s_sex_age_mild_severe.tsv", grp), 
       sep = "\t")

fwrite(res3, file=sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_mild_severe.tsv", grp), 
       sep = "\t")


# look into the genes with NA for pvalues
min_vec = apply(trec1, 1, min)
max_vec = apply(trec1, 1, max)
zero_indexes = which(min_vec == max_vec)
summary(max_vec[zero_indexes])
# they are exactly those with 0 total count for every donor
mean(which(min_vec == max_vec) == which(is.na(pvals0[,1])))



gc()






sessionInfo()
q(save="no")
