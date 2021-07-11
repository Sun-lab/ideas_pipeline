# step11o_DESeq2_ten_p.R
# the main part of the code is carried over from step11g_DESeq2_p.R
# except that this time we do 10 permutations and combine the results
#    - modify the output data frame to hold results from 
#      ten permutations
#    - put results from ten permutations into list

# the main part of the code is carried over from step11b_DESeq2_p.R
# add the permutation part on the diagnosis
# changed the step number and name of the files to save
#    - the labels are permutated after meta_ind is created
#      * the labels in meta is not changed, since it is not used
# add packages and code lines for keeping random seed fixed
library(Matrix)
library(data.table)
library(dplyr)
library(doParallel)
library(doRNG)
library(DESeq2)
#library(svd)

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

print(Sys.getenv("R_HOME"))

RNGkind("L'Ecuyer-CMRG")


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







ori_diagnosis = meta_ind$diagnosis

# set random seed

res0_list = list()
res_list = list()

seed_v = c(5608, 9903, 9968, 7527, 7879, 3760, 2066, 7577, 5926, 9670)

# loop 10 times
for (j in 1:10){

  set.seed(seed_v[j])
  print(paste("j = ", j))
  print(paste("random seed for permutation = ", seed_v[j]))


  meta_ind$diagnosis = sample(ori_diagnosis)



# ------------------------------------------------------------------------
# collect count data
# ------------------------------------------------------------------------


#summary(apply(trec1, 1, median))
#q75 = apply(trec1, 1, quantile, probs=0.75)
#summary(q75)
#table(q75 >= 20)

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

#dd1 = DESeqDataSetFromMatrix(countData = trec1, 
#                             colData = colData,
#                             design = ~ age + sex + Capbatch + PMI + 
#                               RIN + diagnosis)
#dd1 = DESeq(dd1)

#res1 = results(dd1)
#dim(res1)
#res1[1:2,]

  summary(res0)
#summary(res1)

# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

#png(sprintf("figures/step1_DESeq2_%s_compare_pval.png", grp), 
#    width=9, height=3, units="in", res=400)

#par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))
#hist(res0$pvalue, main="without covariates", xlab="p-value")
#hist(res1$pvalue, main="with covariates", xlab="p-value")
#plot(-log10(res0$pvalue), -log10(res1$pvalue), main="-log10(p-value)", 
#     xlab="without covariates", ylab="with covariates", 
#     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
#abline(0, 1, col="darkblue")
#dev.off()

#nms = resultsNames(dd1)
#nms
#nms = nms[-1]

#pvals = matrix(NA, nrow=nrow(trec1), ncol=length(nms))

#for(k in 1:length(nms)){
#  rk = results(dd1, name=nms[k])
#  pvals[,k] = rk$pvalue
#}

#colnames(pvals) = nms
#dim(pvals)
#summary(pvals)

#pdf(sprintf("figures/step1_DESeq2_%s_pval_hist.pdf", grp), width=9, height=9)
#par(mfrow=c(3,3), bty="n", mar=c(5,4,2,1))
#for(k in 1:length(nms)){
#  hist(pvals[,k], main=nms[k], xlab="p-value", breaks=50)
#}

#rbatch = results(dd1, contrast=c("Capbatch","CB6","CB7"))
#dim(rbatch)
#rbatch[1:2,]
#hist(rbatch$pvalue, main="Capbatch_CB7_vs_CB6", xlab="p-value", breaks=50)
#dev.off()

# ------------------------------------------------------------------------
# now we conclude that RMI is not significantly associated with 
# gene expression, and Capbatch CB1 and CB2 are similar, so does 
# Capbatch CB6 and CB7. so we can replace them by sequencing batch
# ------------------------------------------------------------------------

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

#pdf(sprintf("figures/step11o_DESeq2_%s_pval_hist_final_ten_p.pdf", grp), 
#    width=9, height=6)
#par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
#for(k in 1:length(nms)){
#  hist(pvals2[,k], main=nms[k], xlab="p-value", breaks=50)
#}

#plot(-log10(res1$pvalue), -log10(res$pvalue), main="-log10(p-value)", 
#     xlab="with all covariates", ylab="with selected covariates", 
#     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
#abline(0, 1, col="darkblue")

#dev.off()



# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------



#pdf(sprintf("figures/step11o_DESeq2_%s_compare_pval_Seq_ten_p.pdf", grp), 
#    width=9, height=3)
#par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))
#hist(res0$pvalue, main="without covariates", xlab="p-value")
#hist(res$pvalue, main="with covariates (Seq)", xlab="p-value")
#plot(-log10(res0$pvalue), -log10(res$pvalue), main="-log10(p-value)", 
#     xlab="without covariates", ylab="with covariates (Seq)", 
#     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
#abline(0, 1, col="darkblue")
#dev.off()



# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

  dim(res0)
  res0[1:2,]

  dim(res)
  res[1:2,]

  res0 = as.data.frame(res0)
  res  = as.data.frame(res)


  res0_list[[j]] = res0
  res_list[[j]]  = res

}

saveRDS(res0_list, sprintf("res/step11o_DESeq2_%s_no_covariates_ten_p.rds", grp))
saveRDS(res_list, sprintf("res/step11o_DESeq2_%s_adj_covariates_ten_p.rds", grp))

#fwrite(trec1, file=sprintf("res/step11o_DESeq2_%s_edata_ten_p.tsv", grp), 
#            sep = "\t")

#saveRDS(meta, sprintf("res/step11o_DESeq2_%s_meta_cell_ten_p.rds", grp))
#saveRDS(meta_ind, file=sprintf("res/step11o_DESeq2_%s_meta_ind_ten_p.rds", grp))

gc()


sessionInfo()
q(save="no")
