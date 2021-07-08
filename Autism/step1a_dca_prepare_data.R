
library(Matrix)

data.dir = "~/data"

fls = list.files(path=file.path(data.dir, "ct_mtx"), pattern="PFC_")
fls

ctypes = c("L2_3", "L4", "Microglia", "Endothelial", "IN-SV2C", "AST-PP", 
           "IN-VIP", "IN-SST", "IN-PV", "AST-FB", "Oligodendrocytes", 
           "L5_6", "L5_6-CC", "OPC", "Neu-NRGN-II", "Neu-NRGN-I", "Neu-mat")

ct.mtrx = NULL

for(i in 1:length(ctypes)){
  ctp  = ctypes[i]
  
  cat(i, ": ", ctp, date(), "\n")
  
  dat1 = readRDS(file.path(data.dir, "ct_mtx", sprintf("PFC_%s.rds", ctp)))
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

write.csv(ct.mtrx, file.path(data.dir, "PFC_all.csv"))

gc()

sessionInfo()
q(save="no")
