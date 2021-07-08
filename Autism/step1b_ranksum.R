library(Matrix)
library(data.table)
library(dplyr)

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

diagnosis=meta$diagnosis

date()
ranksum_pval=apply(dat1,1,function(x) wilcox.test(x[diagnosis!="Control"],x[diagnosis=="Control"])$p.value)
date()

rm(dat1)

length(ranksum_pval)
ranksum_pval[1:4]

pdf(sprintf("figures/step1b_ranksum_%s_pval_hist_final.pdf", grp), 
    width=3, height=3)
par(mfrow=c(1,1), bty="n", mar=c(5,4,2,1))
hist(ranksum_pval, main="RankSum Test", xlab="p-value", breaks=50)
dev.off()

# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(sprintf("figures/step1b_ranksum_%s_compare_pval_Seq.pdf", grp), 
    width=3, height=3)
par(mfrow=c(1,1), bty="n", mar=c(5,4,1,1))
hist(ranksum_pval, main="RankSum Test", xlab="p-value", breaks=50)
dev.off()
# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

fwrite(data.frame(ranksum_pval), file=sprintf("res/step1b_ranksum_%s.tsv", grp), 
       sep = "\t",row.names=TRUE,col.names=FALSE)

gc()

sessionInfo()
q(save="no")
