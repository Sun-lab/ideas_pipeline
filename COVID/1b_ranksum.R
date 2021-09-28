library(Matrix)
library(data.table)
library(dplyr)

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

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
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
df_donor[1:2, ]
donor2kp = df_donor$Var1[which(df_donor$Freq >= 30)]

meta2kp = meta_no_early[which(meta_no_early$donor %in% donor2kp),]
dim(meta2kp)

cell2kp_index = which(meta$cell %in% meta2kp$cell)


# create a column for cell level label in terms of COVID and control
table(meta2kp$group_per_sample)
meta2kp$diagnosis = meta2kp$group_per_sample
meta2kp$diagnosis[which(meta2kp$group_per_sample == "mild")] = "COVID"
meta2kp$diagnosis[which(meta2kp$group_per_sample == "severe")] = "COVID"
table(meta2kp$diagnosis)
dim(meta2kp)

# select counts of cells to keep
dat1 = dat[, cell2kp_index]
mean(colnames(dat1) == meta2kp$cell)



# ------------------------------------------------------------------------
# run ranksum
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

diagnosis=meta2kp$diagnosis

date()
ranksum_pval=apply(dat1,1,function(x) wilcox.test(x[diagnosis!="control"],x[diagnosis=="control"])$p.value)
date()

rm(dat1)

length(ranksum_pval)
ranksum_pval[1:4]

pdf(sprintf("figures/1b_ranksum_%s_pval_hist_final.pdf", grp), 
    width=3, height=3)
par(mfrow=c(1,1), bty="n", mar=c(5,4,2,1))
hist(ranksum_pval, main="RankSum Test", xlab="p-value", breaks=50)
dev.off()

# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

fwrite(data.frame(ranksum_pval), file=sprintf("res/1b_ranksum_%s.tsv", grp), 
       sep = "\t",row.names=TRUE,col.names=FALSE)

gc()

sessionInfo()
q(save="no")
