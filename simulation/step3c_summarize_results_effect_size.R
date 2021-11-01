
# --------------------------------------------------------------------------
# similar to step3, but summarize results for 360 cells per subject
# and with different effect sizes
# --------------------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(ggpubr)

theme_set(theme_bw())

# --------------------------------------------------------------------------
# check all the result files
# --------------------------------------------------------------------------

res.files = list.files(path="results", 
                       pattern="_10_nctrl_10_ncell_360_\\S+.txt", 
                       full.names=TRUE)
res.files

gg_all = NULL

for(i in 1:length(res.files)){
  file.i = basename(res.files[i])
  ncase  = str_extract(file.i, '(?<=ncase_)\\d+')
  nctrl  = str_extract(file.i, '(?<=nctrl_)\\d+')
  r_mean = str_extract(file.i, '(?<=mean_)(\\d|\\.)+(?=(_|\\.txt))')
  r_var  = str_extract(file.i, '(?<=var_)(\\d|\\.)+(?=(_|\\.txt))')

  config = gsub("pval_", "", file.i)
  config = gsub(".txt", "", config)
  
  pval2 = read.table(res.files[i], header=TRUE, as.is=TRUE)
  dim(pval2)
  pval2[1:2,]

  cal.power <- function(x, geneType){
    tapply(x, geneType, function(x){sum(x < 0.05, na.rm=TRUE)/sum(!is.na(x))})
  }
  
  powers = apply(pval2[,-1], 2, cal.power, geneType=pval2$geneType)
  
  print(config)
  print(powers)
  
  cols2kp = c("ranksum_pval", "mast_pval_glm", "mast_pval_glmer", 
              "deseq2_pval", "PS_zinb_Was", "PS_kde_Was")
  
  powers  = powers[,match(cols2kp, colnames(powers))]
  gg = melt(powers)
  
  names(gg) = c("geneType", "method", "power")
  gg$method = gsub("deseq2_pval", "DEseq2", gg$method)
  gg$method = gsub("PS_zinb_Was", "IDEAS_ZINB",  gg$method)
  gg$method = gsub("PS_kde_Was",  "IDEAS_KDE", gg$method)
  gg$method = gsub("mast_pval_glmer", "MAST_glmer", gg$method)
  gg$method = gsub("mast_pval_glm", "MAST", gg$method)
  gg$method = gsub("ranksum_pval", "Rank-sum", gg$method)
  table(gg$method)
  gg$method = factor(gg$method, 
                     levels = c("Rank-sum", "MAST", "MAST_glmer", 
                                "DEseq2", "IDEAS_ZINB", "IDEAS_KDE"))
  
  g1 = ggplot(subset(gg, geneType %in% c("EE")), 
              aes(x=geneType, y=power, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) + 
    scale_fill_brewer(palette="Paired") + 
    geom_hline(yintercept=0.05, col="red") + 
    ylab("Type I error")
  
  g2 = ggplot(subset(gg, geneType %in% c("meanDE", "varDE")), 
              aes(x=geneType, y=power, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) + 
    scale_fill_brewer(palette="Paired") + 
    geom_hline(yintercept=0.05, col="red")
  
  gg1 = ggarrange(g1, g2, ncol = 2, nrow = 1, widths=c(1.25,2), 
                  common.legend = TRUE, legend = "top")
  
  ggsave(sprintf("figures/power_fold_change_%s.pdf", config), 
         gg1, width=6, height=3)
  
  gg$mean_fold = rep(r_mean, nrow(gg))
  gg$var_fold  = rep(r_var, nrow(gg))
  
  gg_all = rbind(gg_all, gg)
}


# --------------------------------------------------------------------------
# check all the result files
# --------------------------------------------------------------------------


sessionInfo()
q(save = "no")

