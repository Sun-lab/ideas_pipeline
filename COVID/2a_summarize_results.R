
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

# -------------------------------------------------------------------
# read p-values, draw histogram of p-values, and record 
# fisher exact test to compare each pair of method
# -------------------------------------------------------------------

methods = c("DESeq2", "rank_sum", "MAST_glm", "MAST_glmer", 
            "PS_nb_Was", "PS_dca_direct_Was", "PS_saver_direct_Was")

pi0     = matrix(NA, nrow=length(cell_types), ncol=length(methods))
ngenes  = rep(NA, length(cell_types))
or_list = fp_list = pv_list = list()

pdf("figures/2_pval_hist.pdf", width=12, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,3,1))

for(i in 1:length(cell_types)){
  ct1 = cell_types[i]
  p1  = fread(sprintf("res/1e_pvals_%s.tsv", ct1))
  pv_list[[ct1]]= p1
  ngenes[i] = nrow(p1)
  
  odds_ratio = fisher_pvl = matrix(NA, nrow=length(methods), 
                                   ncol=length(methods))
  
  colnames(odds_ratio) = rownames(odds_ratio) = methods
  colnames(fisher_pvl) = rownames(fisher_pvl) = methods
  
  j = 0
  for(m1 in methods){
    j = j + 1
    hist(p1[[m1]], main=paste0(m1, "\n", ct1), xlab="p-value", breaks = 20)
    pi0[i,j] = min(1, 2*mean(p1[[m1]] > 0.5, na.rm=TRUE))
    
    pj = p1[[m1]] < 0.05
    if(sum(pj, na.rm=TRUE) == 0){ next }
    
    for(k1 in methods){
      pk = p1[[k1]] < 0.05
      if(sum(pk, na.rm=TRUE) == 0){ next }
      
      fjk = fisher.test(pj, pk)
      odds_ratio[m1, k1] = fjk$estimate
      fisher_pvl[m1, k1] = fjk$p.value
    }
  }
  plot(0:1,0:1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  or_list[[ct1]] = odds_ratio
  fp_list[[ct1]] = fisher_pvl
  
}

dev.off()


df1 = data.frame(cell_type =  cell_types, n_genes = ngenes)
df1$cell_type = factor(df1$cell_type, levels = df1$cell_type)

g2  = ggplot(data=df1, aes(x=cell_type, y=n_genes)) +
  geom_bar(stat="identity", fill="gray") + 
  ylab("number of genes") + xlab("") + 
  coord_flip() + theme_classic() 

ggsave("figures/2_barplot_n_genes.pdf", g2, width=3, height=3.5)

pdf("figures/2_odds_ratio.pdf", width=6, height=6)
par(mfrow=c(1,1), bty="n", mar=c(5,4,3,1))

for(ct1 in cell_types){
  or1 = or_list[[ct1]]
  or1[which(or1 > 10)] = 10
  g1 = ggcorrplot(or1, tl.srt=90) + ggtitle(ct1) +
    scale_fill_gradient2(limit = c(0,10), low = "blue", 
                         high =  "red", mid = "white", 
                         midpoint = 5)
  print(g1)
}
dev.off()

# -------------------------------------------------------------------
# summarize the proportion of non-nulls
# -------------------------------------------------------------------

colnames(pi0) = methods
rownames(pi0) = cell_types

pi1 = t(1 - pi0)

pi1

rownames(pi1) = gsub("PS_nb_Was", "IDEAS_NB", rownames(pi1))
rownames(pi1) = gsub("PS_dca_direct_Was", "IDEAS_DCA", rownames(pi1))
rownames(pi1) = gsub("PS_saver_direct_Was", "IDEAS_SAVER", rownames(pi1))

summary(c(pi1))

g1 = ggcorrplot(pi1, tl.srt=90)
g1 = g1 + scale_fill_gradient2(limit = c(0,1), low = "blue", 
                               high =  "red", mid = "white", 
                               midpoint = 0.2)

ggsave("figures/2_pi1_complete.pdf", g1, width=5, height=5)

pi1_sub = pi1[-(2:3),, drop = FALSE]
summary(c(pi1_sub))
rownames(pi1_sub) = rownames(pi1_sub)

g1 = ggcorrplot(pi1_sub, tl.srt=90)
g1 = g1 + scale_fill_gradient2(limit = c(0,0.61), low = "blue", 
                               high =  "red", mid = "white", 
                               midpoint = 0.2)

ggsave("figures/2_pi1_PS_Was.pdf", g1, width=4, height=5)

pi1_sub = as.data.frame(round(t(pi1_sub),3))
table(df1$cell_type == rownames(pi1_sub))
pi1_sub
df1 = cbind(df1, pi1_sub)
df1

write.table(df1, file="res/Fig4A_n_genes.txt", sep="\t", quote=FALSE, 
            row.names = FALSE)








gc()

sessionInfo()
q(save="no")
