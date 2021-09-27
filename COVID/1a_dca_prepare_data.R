library(Matrix)

data.dir = "../../ideas_data/COVID/PBMC_10x"

fls = list.files(path=file.path(data.dir, "ct_mtx"))
fls

ctypes = c('Bcells_1', 'Bcells_2', 'Bcells_3','CD163+Monocytes', 
           'CD4+Tcells_1', 'CD4+Tcells_2', 'CD4+Tcells_3', 'CD8+Tcells_1',
           'CD8+Tcells_2', 'CD8+Tcells_3', 'ClassicalMonocytes', 'HLA-DR-S100A+monocytes', 
           'HLA-DR+CD83+Monocytes', 'ImmatureNeutrophils', 'mDCs', 'Megakaryocyte', 
           'mixed', 'Neutrophils', 'NKcells', 'Non-classicalMonocytes', 'pDCs', 
           'Plasmablasts', 'undefined')

ct.mtrx = NULL

for(i in 1:length(ctypes)){
  ctp  = ctypes[i]
  
  cat(i, ": ", ctp, date(), "\n")
  
  dat1 = readRDS(file.path(data.dir, "ct_mtx", sprintf("%s.rds", ctp)))
  dim(dat1)
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

write.csv(ct.mtrx, file.path(data.dir, "all.csv"))

gc()

sessionInfo()
q(save="no")
