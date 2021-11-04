
# --------------------------------------------------------------------------
# similar to step3, but summarize results for 360 cells per subject
# and with different effect sizes
# --------------------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(ggpubr)

theme_set(theme_classic())

# --------------------------------------------------------------------------
# compare results using NB vs. ZINB
# --------------------------------------------------------------------------

res.files = list.files(path="results", 
                       pattern="_10_nctrl_10_ncell_360_\\S+.txt", 
                       full.names=TRUE)
res.files

pval_zinb = read.table(res.files[7], header=TRUE, as.is=TRUE)
dim(pval_zinb)
pval_zinb[1:2,]

pval_nb = read.table(res.files[8], header=TRUE, as.is=TRUE)
dim(pval_nb)
pval_nb[1:2,]

crs = rep(NA, ncol(pval_nb) - 1)

for(k in 2:ncol(pval_nb)){
  crs[k-1] = cor(pval_zinb[[k]], pval_nb[[k]], method="spearman", use="pair")
}
crs

# plot(-log10(pval_zinb[["PS_nb_Was"]]), -log10(pval_nb[["PS_zinb_Was"]]))
# plot(-log10(pval_zinb[["PS_nb_JSD"]]), -log10(pval_nb[["PS_zinb_JSD"]]))

# --------------------------------------------------------------------------
# check all the result files
# --------------------------------------------------------------------------

res.files = res.files[-7]
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

dim(gg_all)
gg_all[1:3,]


p0 = ggplot(subset(gg_all, geneType %in% c("EE")), 
           aes(x=mean_fold, y=power, group=method)) +
  geom_line(aes(color=method)) + ylab("Type I error") + 
  geom_point(aes(color=method, shape=method)) + 
  scale_color_brewer(palette="Dark2")


p0_sub = p0 + ylim(0,0.1)

table(gg_all$method)
gg_sub = subset(gg_all, ! method %in% c("Rank-sum", "MAST", "MAST_glmer"))

p_mean = ggplot(subset(gg_sub, geneType %in% c("meanDE")), 
                aes(x=mean_fold, y=power, group=method)) +
  geom_line(aes(color=method)) + 
  geom_point(aes(color=method, shape=method)) + 
  scale_color_brewer(palette="Dark2")

p_var = ggplot(subset(gg_sub, geneType %in% c("varDE")), 
                aes(x=var_fold, y=power, group=method)) +
  geom_line(aes(color=method)) + 
  geom_point(aes(color=method, shape=method)) + 
  scale_color_brewer(palette="Dark2")

g0 = ggarrange(p0, p0_sub, labels = c("A", "B"), legend = "top", 
               ncol = 2, nrow = 1, common.legend = TRUE)

gp = ggarrange(p_mean, p_var, labels = c("A", "B"), legend = "top", 
          ncol = 2, nrow = 1, common.legend = TRUE)

ggsave(sprintf("figures/effect_size_10_10_360_type_I_error.pdf", config), 
       g0, width=6, height=3)

ggsave(sprintf("figures/effect_size_power_10_10_360.pdf", config), 
       gp, width=6, height=3)

sessionInfo()
q(save = "no")

