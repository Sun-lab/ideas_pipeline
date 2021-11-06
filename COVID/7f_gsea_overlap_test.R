# This file compares the gsea pathways 
# identified by DESeq2, ideas nb, dca_direct,
# mean expression level approach and residual approach


library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

# for now, CD8+Tcells_1 only
cell_types = c("CD8+Tcells_1")

methods = c("DESeq2", "PS_nb_Was", "PS_dca_direct_Was")

# -------------------------------------------------------------------
# read GSEA results
# -------------------------------------------------------------------

f1 <- function(x){
  ww1 = which(x$padj < 0.05)
  if(length(ww1) > 0){
    x1 = x[ww1, c("pathway", "pval", "padj", "NES")]
  }else{
    x1 = NULL
  }
  x1
}



sink("res/7f_GSEA.txt")

ct1 = cell_types[1]
 
cat("\n---------------------------------------------------\n")
cat(ct1)
cat("\n---------------------------------------------------\n")

gsea = readRDS(sprintf("res/5f_gsea_%s.rds", ct1))
gsea = gsea[methods]

gsea_express = readRDS(sprintf("res/5f_gsea_expression_%s.rds",ct1))
gsea['email']  = gsea_express['email']


x1   = lapply(gsea, f1)
x1[sapply(x1, is.null)] = NULL

if(length(x1) > 0){ print(x1) }

sink()

# -------------------------------------------------------------------
# plot GSEA results for all methods
# -------------------------------------------------------------------

ct1 = "CD8+Tcells_1"

all_methods = c(methods, "email")

sign_pathway_list = list()

for (i in 1:length(all_methods)){
  m = all_methods[i]
  cur_method = gsea[[m]]
  cur_padj = cur_method$padj
  cat(paste(m, "  minimal adjusted pvalue: ",
            as.character(min(cur_padj, na.rm = TRUE)), "\n", sep = ""))
  
  row_sign_padj = which(cur_padj < 0.05)
  cat(paste("number of pathways with adjusted pvalue<0.05: ", 
            as.character(length(row_sign_padj)), "\n" ))

  if (length(row_sign_padj)>0){
    cur_method = cur_method[row_sign_padj,]
    sign_pathway_list[[m]] = cur_method$pathway
  }
}


overlap_mat = matrix(NA, ncol = length(all_methods), 
                     nrow = length(all_methods))

for (i in 1:(length(all_methods)-1)){
  for (j in (i+1):length(all_methods)){
    mi = all_methods[i]
    mj = all_methods[j]
    pathways_i = sign_pathway_list[[mi]]
    pathways_j = sign_pathway_list[[mj]]
    len_overlap = length(intersect(pathways_i, pathways_j))
    overlap_mat[i, j] = paste0(len_overlap, "/", length(pathways_i))
    overlap_mat[j, i] = paste0(len_overlap, "/", length(pathways_j)) 
  }
}

first_column = paste0(all_methods, 
                      rep(" (", length(all_methods)),
                      sapply(sign_pathway_list, length), 
                      rep(")", length(all_methods)))

df_overlap = cbind(first_column, overlap_mat)
colnames(df_overlap) = c("methods", all_methods)

write.csv(df_overlap, file = "res/7f_overlap_count_matrix.csv", 
          row.names = FALSE)


or_mat = matrix(NA, ncol = length(all_methods), 
                nrow = length(all_methods))
R_or_mat = matrix(NA, ncol = length(all_methods), 
                  nrow = length(all_methods))
fp_mat = matrix(NA, ncol = length(all_methods), 
                nrow = length(all_methods))

for (i in 1:(length(all_methods)-1)){
  for (j in (i+1):length(all_methods)){
    
    mi = all_methods[i]
    mj = all_methods[j]
    
    cur_padj_i = gsea[[mi]]$padj
    pi = cur_padj_i < 0.05
    if(sum(pi, na.rm=TRUE) == 0){ next }
    
    cur_padj_j = gsea[[mj]]$padj
    pj = cur_padj_j < 0.05
    if(sum(pj, na.rm=TRUE) == 0){ next }
    
    pathways_i = gsea[[mi]]$pathway
    pathways_j = gsea[[mj]]$pathway
    # make the order of the pathways from method j 
    # match those from method i
    cur_padj_j = cur_padj_j[match(pathways_i, pathways_j)]
    pj = cur_padj_j < 0.05
    
    element_2kp = which((!is.na(cur_padj_i))&(!is.na(cur_padj_j)))
    # not sure why, but adding na.rm = TRUE in the following brackets will 
    # at mistakenly introduce an additional 1 if it should sum to 0
    pp = sum((pi & pj)[element_2kp])
    pn = sum((pi & (!pj))[element_2kp])
    np = sum(((!pi) & pj)[element_2kp])
    nn = sum(((!pi) & (!pj))[element_2kp])
    
    fij = fisher.test(pi[element_2kp], pj[element_2kp], greater)
    or_mat[j, i] = 
      format(round((pp*nn)/(pn*np), digits = 2), nsmall = 2)
    R_or_mat[j, i] = format(round(fij$estimate, digits = 2), nsmall = 2)
    fisher_pvalue = fij$p.value
    if (fisher_pvalue >= 0.001){
      fisher_pvalue = 
        format(round(fisher_pvalue, digits = 3), nsmall = 3)
    }else{
      fisher_pvalue = 
        formatC(fisher_pvalue, format = "e", digits = 1)
    }
    fp_mat[j, i] = fisher_pvalue
  }
}



df_or = cbind(first_column, or_mat)
colnames(df_or) = c("methods", all_methods)

write.csv(df_or, file = "res/7f_odds_ratio_matrix.csv", 
          row.names = FALSE)

df_R_or = cbind(first_column, R_or_mat)
colnames(df_R_or) = c("methods", all_methods)

write.csv(df_R_or, file = "res/7f_fisher_odds_ratio_matrix.csv", 
          row.names = FALSE)

df_fp = cbind(first_column, fp_mat)
colnames(df_fp) = c("methods", all_methods)

write.csv(df_fp, file = "res/7f_fisher_pvalue_matrix.csv", 
          row.names = FALSE)






gc()

sessionInfo()
q(save="no")
