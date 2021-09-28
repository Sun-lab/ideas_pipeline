
# ---------------------------------------------------------------
# initial setup
# ---------------------------------------------------------------

library(Matrix)
library(SAVER)
library(pryr)

nCores = 12
# ---------------------------------------------------------------
# load data
# ---------------------------------------------------------------

data.dir = "../../ideas_data/COVID/PBMC_10x"

fls = list.files(path=file.path(data.dir, "ct_mtx"))
fls


ctypes = c('Bcells_1', 'Bcells_2', 'Bcells_3','CD163+Monocytes', 
           'CD4+Tcells_1', 'CD4+Tcells_2', 'CD4+Tcells_3', 'CD8+Tcells_1',
           'CD8+Tcells_2', 'CD8+Tcells_3', 'ClassicalMonocytes', 'HLA-DR-S100A+monocytes', 
           'HLA-DR+CD83+Monocytes', 'ImmatureNeutrophils', 'mDCs', 'Megakaryocyte', 
           'mixed', 'Neutrophils', 'NKcells', 'Non-classicalMonocytes', 'pDCs', 
           'Plasmablasts', 'undefined')

#ctypes = c("Microglia", "Endothelial","Neu-mat")
ct.mtrx = NULL
meta=matrix(ncol=1,nrow=0)

quantile_info_before = matrix(ncol=5,nrow=length(ctypes))
rownames(quantile_info_before) = ctypes
colnames(quantile_info_before) = names(quantile(ct.mtrx))

dim(quantile_info_before)
quantile_info_before[1:2,]

for(i in 1:length(ctypes)){
  ctp  = ctypes[i]
  
  cat(i, ": ", ctp, date(), "\n")
  
  dat1 = readRDS(file.path(data.dir, "ct_mtx", sprintf("%s.rds", ctp)))
  dim(dat1)
  
  quantile_info_before[i,] = quantile(colSums(dat1))
  
  meta=c(meta,rep(ctp,ncol(dat1)))
  
  if(i == 1){
    ct.mtrx = as.matrix(dat1)
  }else{
    ct.mtrx = cbind(ct.mtrx, as.matrix(dat1))
  }
  rm(dat1)
  
  gci = gc()
  print(gc())
  cat("\n")
}

dim(ct.mtrx)
ct.mtrx[1:2,1:5]

table(meta)

# ---------------------------------------------------------------
# saver 
# ---------------------------------------------------------------

dim(ct.mtrx)

# first a test run for memory usage and computational time
date()
mem_used()
gc()
ct.mtrx.saver = saver(ct.mtrx[1:1000, 1:200], 
                      estimates.only = TRUE, ncores=nCores)
date()
mem_used()
gc()


date()
mem_used()
gc()
ct.mtrx.saver = saver(ct.mtrx,estimates.only = TRUE, ncores=nCores)
date()
gc()
mem_used()

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------

quantile_info_after = matrix(ncol=5,nrow=length(ctypes))
rownames(quantile_info_after) = ctypes
colnames(quantile_info_after) = names(quantile(ct.mtrx.saver))

for(i in 1:length(ctypes)){
  ctp  = ctypes[i]
  
  cat(i, ": ", ctp, date(), "\n")
  
  dat1=ct.mtrx.saver[,which(meta==ctp)]
   
  quantile_info_after[i,] = quantile(colSums(dat1)) 
  print(dim(dat1))
  saveRDS(dat1,file.path(data.dir, "ct_mtx", sprintf("saver_%s.rds", ctp)))
  
  gci = gc()
  print(gc())
  cat("\n")
}

dat1[1:3,1:2]
quantile_info_before
quantile_info_after

sessionInfo()
mem_used()
gc()
q(save="no")
