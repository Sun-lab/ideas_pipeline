library(data.table)
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(ggrepel)
library(stringr)
library(dplyr)
library(Matrix)
library(segmented)
library(Seurat)
library(maftools)
theme_set(theme_bw())

# use the data from cohort 1 PBMC 10x 

path = "../../ideas_data/COVID/"
folder_2 = "dataset-952687f71ef34322a850553c4a24e82e/"
fnm_2  = "seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"
PBMC_cohort1_10x_jonas_FG  = readRDS(paste0(path, folder_2, fnm_2))

# data exploration

PBMC_cohort1_10x_jonas_FG 

print(slotNames(PBMC_cohort1_10x_jonas_FG))

length(slotNames(PBMC_cohort1_10x_jonas_FG))

dat_2 = PBMC_cohort1_10x_jonas_FG@assays$RNA
dat_2

print(dim(dat_2))
dat_2[1:6, 1:2]

cell_info = PBMC_cohort1_10x_jonas_FG@meta.data
print(dim(cell_info))
print(names(cell_info))

cell_info[1:2, 1:6]

table(cell_info$id.celltype)

sort(table(cell_info$id.celltype))

table(cell_info$sampleID, cell_info$id.celltype)

table(cell_info$disease_stage)

table(cell_info$disease_stage, cell_info$id.celltype)

table(cell_info$donor, cell_info$disease_stage)

id = which(cell_info$id.celltype == '12: CD8+ T cells_1')

table(cell_info[id,]$donor, cell_info[id,]$disease_stage)

cell_info[1:2, 1:10]

cell_info[1:2, 11:20]

cell_info[1:2, 21:30]

cell_info[1:2, 31:39]

table(cell_info[id,]$donor, cell_info[id,]$group_per_sample)

table(cell_info$group_per_sample)

df2_meta_ind = distinct(cell_info[, c("donor", "group_per_sample")])

df2_meta_ind

length(unique(df2_meta_ind$donor))

table(df2_meta_ind$group_per_sample)

# Recover the count matrix from assays RNA matrix

print(dim(dat_2))

mean(colnames(dat_2) == rownames(cell_info))

# this step is not necessary, since it has been verified that the order of col names in 
# assay RNA matrix is the same as that in the row names for cell info
meta_cell = cell_info[match(colnames(dat_2), rownames(cell_info)),]
meta_cell[1:6, 1:2]

mean(cell_info$nCount_RNA == meta_cell$nCount_RNA)

cell_type_list =  unique(c(cell_info$id.celltype))
print(length(cell_type_list))
cell_type_list

cell_type = as.character(cell_type_list[1])

length(which(cell_info$id.celltype == cell_type))

table(cell_info$id.celltype)

# Write out the meta data file to tsv

write.table(PBMC_cohort1_10x_jonas_FG@meta.data, file = "../../ideas_data/COVID/PBMC_10x/meta.tsv", row.names=TRUE,  sep="\t")

# Write out recovered count matrix for each cell type

data_dir = "../../ideas_data/COVID/PBMC_10x/"

for (i in c(1:length(cell_type_list))){

    cell_type = as.character(cell_type_list[i])
    cell_ids = which(cell_info$id.celltype == cell_type)

    cur_matrix = dat_2[, cell_ids]
    cur_meta = cell_info[cell_ids, ]

    cur_exp_matrix = exp(cur_matrix)
    cur_exp_minus_1_matrix = cur_exp_matrix - 1
    cur_exp_minus_1_div10000_matrix = cur_exp_minus_1_matrix/10000
    
    cur_recover = matrix(nrow = nrow(cur_matrix), ncol = ncol(cur_matrix))
    
    UMI = cur_meta$nCount_RNA
    cur_recover = t(t(cur_exp_minus_1_div10000_matrix)*UMI)
    
    cur_recover_sparse = Matrix(cur_recover, sparse = TRUE) 
    colnames(cur_recover_sparse) = colnames(cur_matrix)
    rownames(cur_recover_sparse) = rownames(cur_matrix)
    cell_type_parse = strsplit(cell_type, " ")[[1]]
    cell_type_write = paste(cell_type_parse[-1], collapse = "")
    saveRDS(cur_recover_sparse, file = paste(data_dir, "ct_mtx/", cell_type_write, ".rds", sep = ""))
    
    print(paste("Done with cell type: ", cell_type, sep = ""))

}


