# the current version of DCA (as of Sept.28, 2021) does not 
# provide the mean_norm.tsv as one of the outputs by default
# this code recovers mean_norm.tsv making use of 
# (1) original count matrix to compute size factors
# (2) mean.tsv



library(Matrix)
library(data.table)

dca.dir = "../../ideas_data/COVID/PBMC_10x/dca_zinb_all"
data.dir = "../../ideas_data/COVID/PBMC_10x"

# Load the original count matrix used for DCA training

all_counts = read.csv(file.path(data.dir, "all.csv"), 
                      header = TRUE, sep = ",", check.names=FALSE)
dim(all_counts)

all_counts[1:6, 1:3]

col_sums = apply(all_counts[, 2:ncol(all_counts)], 2, sum)
length(col_sums)

as.numeric(col_sums[1:10])

med_rd = median(col_sums)
med_rd

size_factors = as.numeric(col_sums)/med_rd
size_factors[1:10]

summary(col_sums)

median(size_factors)



# load mean matrix

mean_mat = fread(file.path(dca.dir, "mean.tsv"), 
                 header = TRUE, sep = "\t")
dim(mean_mat)

mean(all_counts[, 1] == unlist(mean_mat[, 1]))

# Recover mean norm matrix.

mean_norm_mat = t(t(mean_mat[, 2:ncol(mean_mat)])/size_factors)
dim(mean_norm_mat)

mean_names = names(mean_mat)
mean_names[1:10]

mean_norm_names = colnames(mean_norm_mat)
mean_norm_names[1:9]

mean_mat[1:6, 2:3]
mean_norm_mat[1:6, 1:2]

mean(mean_names[2:ncol(mean_mat)] == mean_norm_names)

# test recovery for mean_norm_mat, using certain cells
j = 10000
summary(mean_mat[, (j+1):(j+1)]/as.numeric(unlist(mean_norm_mat[, j:j])))
sum(as.numeric(unlist(mean_norm_mat[, j:j])))
sum(as.numeric(unlist(mean_norm_mat[, j:j]))) * sum(all_counts[, (j+1)])
sum(as.numeric(unlist(mean_mat[, (j+1):(j+1)]))) * med_rd

rnames = rownames(mean_mat)
length(rnames)

mean(colnames(all_counts)[2:ncol(all_counts)] == names(mean_mat)[2:ncol(mean_mat)])


mean_norm_dt = as.data.table(mean_norm_mat)
dim(mean_norm_dt)
mean_norm_dt[1:6, 1:3]

mean_norm_dt[, V1:= mean_mat[, 1]]
dim(mean_norm_dt)

mean_norm_dt[1:6, 99048:99050]
mean_norm_mat[1:6, 99047:99049]

setcolorder(mean_norm_dt, "V1")

mean_norm_dt[1:6, 1:3]
mean_norm_mat[1:6, 1:2]

fwrite(mean_norm_dt, file = file.path(dca.dir, "mean_norm.tsv"), sep="\t")




