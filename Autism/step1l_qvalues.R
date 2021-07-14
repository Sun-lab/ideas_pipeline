# convert pvals to qvals
# count # of significant genes from qvalue based on thresholds
# threshold = c(0.1, 0.2, 0.3, 0.4)

# should be right after step1g


# ========================================================================
# libraries and path
# ========================================================================

library(MASS)
library(data.table)
library(doParallel)
library(doRNG)
library(qvalue)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(fgsea)
library(stringr)



ctypes = c("L2_3", "L4", "Microglia", "Endothelial", "IN-SV2C", "AST-PP", 
           "IN-VIP", "IN-SST", "IN-PV", "AST-FB", "Oligodendrocytes", 
           "L5_6", "L5_6-CC", "OPC", "Neu-NRGN-II", "Neu-NRGN-I", "Neu-mat")
ctypes = sort(ctypes)
length(ctypes)

for (threshold in c(0.1, 0.2, 0.3, 0.4)){
  # the following line contains a hard coded ncol = 10
  # may change to make it more flexible
  sig_count_mat = matrix(ncol = 17, nrow = length(ctypes))
  
  check_length <- function(v){
    length(which(v < threshold))
  }
  
  for (type_ind in 1:length(ctypes)){
    
    grp = paste("PFC_", ctypes[type_ind], sep = "")
    print(sprintf("grp = %s", grp))
    pvals = fread(sprintf("res/step1e_pvals_%s.tsv", grp))
    
    methods = names(pvals)[-1]
    dim(pvals)
    head(pvals)
    
    qvals = matrix(ncol = ncol(pvals), nrow = nrow(pvals))
    qvals = as.data.frame(qvals)
    names(qvals) = names(pvals)
    qvals$gene   = pvals$gene
    
    pi_list = list()
    
    for(method in methods){
      pval_v = pvals[[method]]
      pi_list[[method]] = min(1, 2*sum(pval_v > 0.5, na.rm = TRUE)/sum(!is.na(pval_v)))
      if(pi_list[[method]] < 0.5){
        pi_list[method] = 0.5
      }
      qval   = qvalue(pval_v, pi0=pi_list[[method]])
      qvals[[method]] = qval$qvalues
    }
    
    sig_count_mat[type_ind, 2:ncol(sig_count_mat)] = apply(qvals[, -1], 2, check_length)
    
    if (threshold == 0.4){
      file.name.q = sprintf("res/step1l_qvals_%s.tsv", grp)
      fwrite(qvals, file = file.name.q, sep = "\t")
    }  
  }
  
  sig_count_df = as.data.frame(sig_count_mat)
  names(sig_count_df) = c("cell_type", methods)
  sig_count_df$cell_type = ctypes
  
  file.name = sprintf("res/step1l_counts_sig_qval_original_ctypes_threshold_%s.tsv",
                      as.character(threshold))
  fwrite(sig_count_df, file = file.name, sep = "\t")
  
}




sessionInfo()
q(save="no")
