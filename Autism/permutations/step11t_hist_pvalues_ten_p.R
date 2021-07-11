# this file plots the histograms of pvalues from original run 
# and those from the ten permutations

# main part of the code is carried over from
# step11m_hist_pvalues_p.R

# this file plots the histograms of pvalues from original run 
# and those from the permutated run

library(MASS)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)



ctypes = c("L2_3", "L4", "Microglia", "Endothelial", "IN-SV2C", "AST-PP", 
           "IN-VIP", "IN-SST", "IN-PV", "AST-FB", "Oligodendrocytes", 
           "L5_6", "L5_6-CC", "OPC", "Neu-NRGN-II", "Neu-NRGN-I", "Neu-mat")

# don't know why, but the for loop does not work for 
# writing the pdf files
# have to run it step by step
i = 0
#for (i in 1:length(ctypes)){
  i = i + 1
  grp = paste("PFC_", ctypes[i], sep = "")

  pvals   = fread(sprintf("res/step11e_pvals_%s.tsv", grp))
  pvals_p = fread(sprintf("res/step11r_pvals_%s_ten_p.tsv", grp))

  gh = list()
  
  methods = names(pvals)[2:10]

  for(nm in methods){
    nm_left = paste(nm, ", original", sep = "")
    gh[[nm_left]] = ggplot(pvals, aes_string(x = nm)) + 
      labs(title = nm_left) + 
      geom_histogram(color = "darkblue", fill = "lightblue", 
                     breaks = seq(0,1,by = 0.02))
    
    nm_right = paste(nm, ", ten permutations", sep = "")
    #df = data.frame(x = pvals[[nm]], y = qvals[[nm]])
    gh[[nm_right]] = ggplot(pvals_p, aes_string(x = nm)) + 
      labs(title = nm_right) + 
      geom_histogram(color = "darkblue", fill = "lightblue", 
                     breaks = seq(0,1,by = 0.02))
    }

   pdf(sprintf("figures/step11t_ori_pval_vs_ten_permutated_hist_ori_ctypes_%s.pdf", grp), 
    width=12, height=10)
   ggarrange(plotlist=gh, ncol = 4, nrow = 5)
   dev.off()

#}