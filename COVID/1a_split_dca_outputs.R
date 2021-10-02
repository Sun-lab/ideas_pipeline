
library(data.table)
library(ggplot2)
library(emdbook)
library(stringr)



data.dir = "../../ideas_data/COVID/PBMC_10x"


# ---------------------------------------------------------------------
# a function to split the statistics from DCA into individual files
# for each cell type separately
# ---------------------------------------------------------------------

split_cluster <- function(dat, cell_info, statistics){
  cell_ids = names(dat)[-1]
  length(cell_ids)
  cell_ids[1:5]
  
  # the current version of dca (as of Sept.26, 2021) does not provide
  # meaningful headers (cell names) for the dropout matrix
  # so cannot verify the matching between dropout matrix column and
  # cells in meta info for now
  #stopifnot(all(cell_ids %in% cell_info$cell))
  
  celltypes = cell_info$id.celltype[match(cell_ids, cell_info$cell)]
  table(celltypes)
  
  ucelltypes = unique(celltypes)
  
  for(c1 in ucelltypes){
    # index +1 since dca_mean has one extra column of gene id
    c1 = as.character(c1)   
    w2kp = which(celltypes == c1) + 1
    dat_celltype = dat[,c(1, w2kp),with=FALSE]
    c1 = strsplit(c1, " ")[[1]]
    c1 = paste(c1[-1], collapse = "")
    fnm = file.path(data.dir, 
                    sprintf("dca_zinb_celltypes/%s_%s.tsv", c1, statistics))
    fwrite(dat_celltype, file = fnm, sep = "\t")
    # also tried to save to rds, it seems the file size is not much difference
  }
}

# ---------------------------------------------------------------------
# read in cell information and read count
# ---------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"), stringsAsFactors=TRUE)
dim(cell_info)
cell_info[1:2,]

table(cell_info$id.celltype)


# ---------------------------------------------------------------------
# read in dca output and split data into different cell types
# ---------------------------------------------------------------------


# first is the intermediate level output mean_norm.tsv, before scaling
# back to the raw scale
dca_mean_norm  = fread(file.path(data.dir, "dca_zinb_all/mean_norm.tsv"))
gc()
dim(dca_mean_norm)
dca_mean_norm[1:2,1:5]
# can only do the following verification directly for mean_norm
stopifnot(all(names(dca_mean_norm)[-1] %in% cell_info$cell))

split_cluster(dca_mean_norm, cell_info, "mean_norm")

gc()


# second is the pi parameter
dca_pi = fread(file.path(data.dir, "dca_zinb_all/dropout.tsv"))
gc()
dim(dca_pi)
dca_pi[1:2,1:5]
# set column names of dropout matrix 
# to follow those of mean_norm matrix
names(dca_pi) = names(dca_mean_norm)
  
split_cluster(dca_pi, cell_info, "pi")

gc()


# third is the dispersion parameter
dca_theta = fread(file.path(data.dir, "dca_zinb_all/dispersion.tsv"))
gc()
dim(dca_theta)
dca_theta[1:2,1:5]
# set column names of dispersion matrix 
# to follow those of mean_norm matrix
names(dca_theta) = names(dca_mean_norm)

split_cluster(dca_theta, cell_info, "dispersion")

gc()







sessionInfo()
q(save="no")


