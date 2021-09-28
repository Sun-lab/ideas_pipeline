

data.dir = "../../ideas_data/COVID/PBMC_10x"


library(data.table)
library(ggplot2)
library(emdbook)
library(stringr)

# ---------------------------------------------------------------------
# a function to split the statistics from DCA into individual files
# for each cell type separately
# ---------------------------------------------------------------------

split_cluster <- function(dat, cell_info, statistics){
  cell_ids = names(dat)[-1]
  length(cell_ids)
  cell_ids[1:5]
  
  stopifnot(all(cell_ids %in% cell_info$cell))
  
  clusters = cell_info$id.celltype[match(cell_ids, cell_info$cell)]
  table(clusters)
  
  uClusters = unique(clusters)
  
  for(c1 in uClusters){
    # index +1 since dca_mean has one extra column of gene id
    w2kp = which(clusters == c1) + 1
    dat_clust = dat[,c(1, w2kp),with=FALSE]
    cell_type_parse = strsplit(c1, " ")[[1]]
    cell_type_write = paste(cell_type_parse[-1], collapse = "")
    fnm = file.path(data.dir, sprintf("dca/%s_%s.tsv",  cell_type_write, statistics))
    fwrite(dat_clust, file = fnm, sep = "\t")
    # also tried to save to rds, it seems the file size is not much difference
  }
}

# ---------------------------------------------------------------------
# read in cell information and read count
# ---------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"), stringsAsFactors=TRUE)
dim(cell_info)
cell_info[1:2,]

table(cell_info$cluster)
table(cell_info$cluster[cell_info$region=="PFC"])


# ---------------------------------------------------------------------
# read in dca output and split data into different cell types
# ---------------------------------------------------------------------

# first is the intermediate level output mean_norm.tsv, before scaling
# back to the raw scale
dca_mean_norm  = fread(file.path(data.dir, "dca_zinb_all/mean_norm.tsv"))
gc()
dim(dca_mean_norm)
dca_mean_norm[1:2,1:5]

split_cluster(dca_mean_norm, cell_info, "mean_norm")

gc()

# second is the pi parameter
dca_pi = fread(file.path(data.dir, "dca_zinb_all/pi.tsv"))
gc()
dim(dca_pi)
dca_pi[1:2,1:5]

split_cluster(dca_pi, cell_info, "pi")

gc()

# third is the dispersion parameter
dca_theta = fread(file.path(data.dir, "dca_zinb_all/dispersion.tsv"))
gc()
dim(dca_theta)
dca_theta[1:2,1:5]

split_cluster(dca_theta, cell_info, "dispersion")

gc()







sessionInfo()
q(save="no")


