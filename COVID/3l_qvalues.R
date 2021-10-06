# compared to 1l, the results delt with in this version consider more genes
# by setting filtering criterion to keep genes
# appearning in at last 10% of the cells


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


# for now, CD8+Tcells_1 only
ctypes = c("CD8+Tcells_1")
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
    
    grp = ctypes[type_ind]
    print(sprintf("grp = %s", grp))
    pvals = fread(sprintf("res/3e_pvals_%s.tsv", grp))
    
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
      file.name.q = sprintf("res/3l_qvals_%s.tsv", grp)
      fwrite(qvals, file = file.name.q, sep = "\t")
    }  
  }
  
  sig_count_df = as.data.frame(sig_count_mat)
  names(sig_count_df) = c("cell_type", methods)
  sig_count_df$cell_type = ctypes
  
  file.name = sprintf("res/3l_counts_sig_qval_original_ctypes_threshold_%s.tsv",
                      as.character(threshold))
  fwrite(sig_count_df, file = file.name, sep = "\t")
  
}




sessionInfo()
q(save="no")
