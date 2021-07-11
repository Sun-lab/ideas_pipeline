
library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())
data.dir  = "./data"

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

# -------------------------------------------------------------------
# load MAST results after permutations
# -------------------------------------------------------------------

grp = "PFC_L2_3"
mast_glm = fread(sprintf("res/step1b_MAST_perm_%s_glm.tsv", grp))

dim(mast_glm)
mast_glm[1:2,]
summary(mast_glm$V2)


full_dat1  = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
full_genes = rownames(full_dat1)  

n.zeros = rowSums(full_dat1 == 0)
summary(n.zeros)
w2kp = which(n.zeros < 0.8*ncol(full_dat1))
length(w2kp)
summary(mast_glm$V2[w2kp])


gc()

sessionInfo()
q(save="no")
