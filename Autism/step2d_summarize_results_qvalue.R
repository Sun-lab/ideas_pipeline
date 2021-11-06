# this code compares genes thresholded by qvalue
# from each method

library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

# -------------------------------------------------------------------
# read p-values, draw histogram of p-values, and record 
# fisher exact test to compare each pair of method
# -------------------------------------------------------------------

methods = c("DESeq2", "rank_sum", "MAST_glm", "MAST_glmer", 
            "PS_nb_Was", "PS_dca_direct_Was", "PS_saver_direct_Was")

or_list = fp_list = list()
n_overlap_list = prop_overlap_list = list()


for(i in 1:length(cell_types)){
  ct1 = cell_types[i]
  q1  = fread(sprintf("res/step1l_qvals_%s.tsv", ct1))

  odds_ratio = fisher_pvl = matrix(NA, nrow=length(methods), 
                                   ncol=length(methods))
  n_overlap = prop_overlap = matrix(NA, nrow=length(methods), 
                                    ncol=length(methods))
  
  colnames(odds_ratio) = rownames(odds_ratio) = methods
  colnames(fisher_pvl) = rownames(fisher_pvl) = methods
  
  colnames(n_overlap) = rownames(n_overlap) = methods
  colnames(prop_overlap) = rownames(prop_overlap) = methods
  
  j = 0
  for(m1 in methods){
    j = j + 1
    # use q value cutoff 0.1 instead of pvalue cutoff 0.05
    qj = q1[[m1]] < 0.1
    if(sum(qj, na.rm=TRUE) == 0){ next }
    
    for(k1 in methods){
      # use q value cutoff 0.1 instead of pvalue cutoff 0.05
      qk = q1[[k1]] < 0.1
      if(sum(qk, na.rm=TRUE) == 0){ next }
      # here the alternative is modified to greater
      fjk = fisher.test(qj, qk, alternative = "greater")
      odds_ratio[m1, k1] = fjk$estimate
      fisher_pvl[m1, k1] = fjk$p.value
      n_overlap[m1, k1] = sum(qj & qk, na.rm = TRUE)
      prop_overlap[m1, k1] = sum(qj & qk, na.rm = TRUE)/sum(qj, na.rm = TRUE)
    }
  }
  #plot(0:1,0:1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  or_list[[ct1]] = odds_ratio
  fp_list[[ct1]] = fisher_pvl
  n_overlap_list[[ct1]] = n_overlap
  prop_overlap_list[[ct1]] = round(prop_overlap, digits = 3)
  
}



pdf("figures/step2d_odds_ratio_qvalue_threshold.pdf", width=6, height=6)
par(mfrow=c(1,1), bty="n", mar=c(5,4,3,1))

for(ct1 in cell_types){
  or1 = or_list[[ct1]]
  or1[which(or1 > 10)] = 10
  g1 = ggcorrplot(or1, tl.srt=90) + ggtitle(gsub("PFC_", "", ct1)) +
    scale_fill_gradient2(limit = c(0,10), low = "blue", 
                         high =  "red", mid = "white", 
                         midpoint = 1)
  print(g1)
}
dev.off()



pdf("figures/step2d_fisher_exact_neglogpvals_qvalue_threshold.pdf", width=6, height=6)
par(mfrow=c(1,1), bty="n", mar=c(5,4,3,1))

for(ct1 in cell_types){
  fp1 = fp_list[[ct1]]
  neg_log_fp1 = -log10(fp1)
  neg_log_fp1[which(neg_log_fp1 > 10)] = 10
  g1 = ggcorrplot(neg_log_fp1, tl.srt=90) + ggtitle(gsub("PFC_", "", ct1)) +
    scale_fill_gradient2(limit = c(0,10), low = "blue", 
                         high =  "red", mid = "white", 
                         midpoint = 2)
  print(g1)
}
dev.off()


for(ct1 in cell_types){
  filename_n_overlap = sprintf("res/step2d_n_overlap_%s.csv", ct1)
  write.csv(n_overlap_list[[ct1]], file = filename_n_overlap, 
            row.names = TRUE)
  filename_prop_overlap = sprintf("res/step2d_prop_overlap_%s.csv", ct1)
  write.csv(prop_overlap_list[[ct1]], file = filename_prop_overlap,
            row.names = TRUE)  
}







gc()

sessionInfo()
q(save="no")
