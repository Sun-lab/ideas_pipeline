
# --------------------------------------------------------------------------
# similar to step3, but with results from all IDEAS options
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

res.files = list.files(path="results", pattern="pval_ncase_\\S+.txt", 
                       full.names=TRUE)
res.files

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
  
  cols2rm = c("mast_pval_glm")
  powers  = powers[,which(!colnames(powers) %in% cols2rm)]
  gg = melt(powers)
  
  names(gg) = c("geneType", "method", "power")
  gg$method = gsub("deseq2_pval", "DEseq2", gg$method)
  gg$method = gsub("mast_pval_glmer", "MAST", gg$method)
  
  gg$method = factor(gg$method, levels = c("DEseq2", "MAST", 
    "KR_kde_JSD", "PS_kde_JSD", "KR_kde_Was", "PS_kde_Was", 
    "KR_zinb_JSD", "PS_zinb_JSD", "KR_zinb_Was", "PS_zinb_Was"))
  
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
  
  ggsave(sprintf("figures/power_IDEAS_all_%s.pdf", config), 
         gg1, width=6, height=3)
}

sessionInfo()
q(save = "no")

