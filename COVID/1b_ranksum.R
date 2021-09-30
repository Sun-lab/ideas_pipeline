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
# filter out cells from control samples
# ------------------------------------------------------------------------

table(meta$donor, meta$disease_stage)

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

diagnosis=meta2kp$group_per_sample

date()
ranksum_pval=apply(dat1,1,function(x) wilcox.test(x[diagnosis!="severe"],x[diagnosis=="severe"])$p.value)
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
