# for each cell type, this code plots the negative log of 
# fisher exact test pvalues between any two methods


library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())

# -------------------------------------------------------------------
# read in cell type information
# -------------------------------------------------------------------

cell_types = c("CD8+Tcells_1")

# -------------------------------------------------------------------
# read p-values, draw histogram of p-values, and record 
# fisher exact test to compare each pair of method
# -------------------------------------------------------------------

methods = c("DESeq2", "rank_sum", "MAST_glm", "MAST_glmer", 
            "PS_nb_Was", "PS_dca_direct_Was", "PS_saver_direct_Was")

pi0     = matrix(NA, nrow=length(cell_types), ncol=length(methods))
ngenes  = rep(NA, length(cell_types))
or_list = fp_list  = list()


for(i in 1:length(cell_types)){
  ct1 = cell_types[i]
  p1  = fread(sprintf("res/1e_pvals_%s.tsv", ct1))
  
  odds_ratio = fisher_pvl = matrix(NA, nrow=length(methods), 
                                   ncol=length(methods))
  
  colnames(odds_ratio) = rownames(odds_ratio) = methods
  colnames(fisher_pvl) = rownames(fisher_pvl) = methods
  
  j = 0
  for(m1 in methods){
    j = j + 1
    pi0[i,j] = min(1, 2*mean(p1[[m1]] > 0.5, na.rm=TRUE))
    
    pj = p1[[m1]] < 0.05
    if(sum(pj, na.rm=TRUE) == 0){ next }
    
    for(k1 in methods){
      pk = p1[[k1]] < 0.05
      if(sum(pk, na.rm=TRUE) == 0){ next }
      # here the alternative is modified to greater
      fjk = fisher.test(pj, pk, alternative = "greater")
      odds_ratio[m1, k1] = fjk$estimate
      fisher_pvl[m1, k1] = fjk$p.value
    }
  }
  plot(0:1,0:1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  or_list[[ct1]] = odds_ratio
  fp_list[[ct1]] = fisher_pvl
  
}




pdf("figures/2c_fisher_exact_neg_log_pvals.pdf", width=6, height=6)
par(mfrow=c(1,1), bty="n", mar=c(5,4,3,1))

for(ct1 in cell_types){
  fp1 = fp_list[[ct1]]
  neg_log_fp1 = -log(fp1)
  neg_log_fp1[which(neg_log_fp1 > 10)] = 10
  g1 = ggcorrplot(neg_log_fp1, tl.srt=90) + ggtitle(gsub("PFC_", "", ct1)) +
    scale_fill_gradient2(limit = c(0,10), low = "blue", 
                         high =  "red", mid = "white", 
                         midpoint = 2)
  print(g1)
}
dev.off()





gc()

sessionInfo()
q(save="no")
