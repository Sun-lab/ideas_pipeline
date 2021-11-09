
library(Matrix)
library(data.table)
library(dplyr)
library(DESeq2)


data.dir = "../../ideas_data/COVID/PBMC_10x"
label.dir = "../../ideas_data/COVID"
  
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
# read in and attach perm information
# ------------------------------------------------------------------------

perm_table = fread(file.path(label.dir, "1b_permlabel.tsv"),header=TRUE)
dim(perm_table)
perm_table[1:2,1:3]

# check that each individual has one label
table(tapply(perm_table$p1, meta2kp$donor, 
             function(v){length(unique(v))}))



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


# get label after permutation
iperm = 1
cur_perm_cell_label = data.frame(perm_table)[, iperm]
cur_perm_ind_label = rep(NA, nrow(meta_ind))

for (i in 1:nrow(meta_ind)){
  wi = which(meta2kp$donor == meta_ind$donor[i])
  cur_perm_ind_label[i] = cur_perm_cell_label[wi[1]]
}


meta_ind[, 2] = cur_perm_ind_label

# only keep the 14 individuals with mild or severe label
# after permutation
ind2kp = which(cur_perm_ind_label != "")


colData = meta_ind[ind2kp,]
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

for (i in ind2kp){
  wi = which(meta2kp$donor == meta_ind$donor[i])
  total_rd = c(total_rd, sum(apply(dat1[,wi], 2, sum)))
}

colData$totalrd = total_rd



# fourth, include both donor level total read depth, sex and age

dd3 = DESeqDataSetFromMatrix(countData = trec1[, ind2kp], 
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




pdf(sprintf("figures/1b_DESeq2_%s_pval_hists_perm.pdf", grp), 
    width=4.5, height=3.5)
k = length(nm3)
hist(pvals3[,k], main="log(totalrd)+sex+age", xlab="p-value", breaks=20)
dev.off()


pdf(sprintf("figures/1b_DESeq2_%s_pval_hists_full_perm.pdf", grp), 
    width=9, height=7)
par(mfrow=c(2,2), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nm3)){
  hist(pvals3[,k], main=nm3[k], xlab="p-value", breaks=50)
}

dev.off()



# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

dim(res3)
res3[1:2,]

res3 = as.data.frame(res3)


fwrite(res3, file=sprintf("res/1b_DESeq2_%s_logtotalrd_sex_age_perm.tsv", grp), 
       sep = "\t")





gc()






sessionInfo()
q(save="no")
