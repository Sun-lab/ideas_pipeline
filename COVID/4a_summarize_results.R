# compared with 2a, this version deals with the results when 
# keep all genes appearing in at least 10% of cells



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

pdf("figures/4_pval_hist.pdf", width=12, height=6)
par(mfrow=c(2,4), bty="n", mar=c(5,4,3,1))

for(i in 1:length(cell_types)){
  ct1 = cell_types[i]
  p1  = fread(sprintf("res/3e_pvals_%s.tsv", ct1))
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

ggsave("figures/4_barplot_n_genes.pdf", g2, width=3, height=3.5)

pdf("figures/4_odds_ratio.pdf", width=6, height=6)
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

ggsave("figures/4_pi1_complete.pdf", g1, width=5, height=5)

pi1_sub = pi1[-(2:3),, drop = FALSE]
summary(c(pi1_sub))
rownames(pi1_sub) = rownames(pi1_sub)

g1 = ggcorrplot(pi1_sub, tl.srt=90)
g1 = g1 + scale_fill_gradient2(limit = c(0,1), low = "blue", 
                               high =  "red", mid = "white", 
                               midpoint = 0.2)

ggsave("figures/4_pi1_PS_Was.pdf", g1, width=4, height=5)

pi1_sub = as.data.frame(round(t(pi1_sub),3))
table(df1$cell_type == rownames(pi1_sub))
pi1_sub
df1 = cbind(df1, pi1_sub)
df1

write.table(df1, file="res/4_Fig4A_n_genes.txt", sep="\t", quote=FALSE, 
            row.names = FALSE)


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


gesaL = list()

sink("res/4_GSEA.txt")

for(ct1 in cell_types){
  
  cat("\n---------------------------------------------------\n")
  cat(ct1)
  cat("\n---------------------------------------------------\n")
  
  gsea = readRDS(sprintf("res/3f_gsea_%s.rds", ct1))
  gsea = gsea[methods]
  
  x1   = lapply(gsea, f1)
  x1[sapply(x1, is.null)] = NULL
  
  if(length(x1) > 0){ print(x1) }
  
  gesaL[[ct1]] = gsea
}
sink()



# -------------------------------------------------------------------
# check GSEA results for one cell type
# -------------------------------------------------------------------

plot_gsea <- function(gseaL, ct1){
  # ct1 = "CD8+Tcells_1"
  DESeq2    = gesaL[[ct1]]$DESeq2
  IDEAS_NB  = gesaL[[ct1]]$PS_nb_Was
  IDEAS_DCA = gesaL[[ct1]]$PS_dca_direct_Was
  
  dim(DESeq2)
  dim(IDEAS_NB)
  dim(IDEAS_DCA)
  
  pathways1 = intersect(DESeq2$pathway, IDEAS_NB$pathway)
  pathways2 = intersect(DESeq2$pathway, IDEAS_DCA$pathway)
  pathways = sort(union(pathways1, pathways2))
  length(pathways)
  
  mat1 = match(pathways, DESeq2$pathway)
  mat2 = match(pathways, IDEAS_NB$pathway)
  mat3 = match(pathways, IDEAS_DCA$pathway)
  
  gd = data.frame(pathway = pathways, 
                  DESeq2 = -log10(DESeq2$pval[mat1]),
                  IDEAS_NB = -log10(IDEAS_NB$pval[mat2]),
                  IDEAS_DCA = -log10(IDEAS_DCA$pval[mat3]))
  
  g1 = ggplot(gd, aes(x=DESeq2, y=IDEAS_NB)) + geom_point(alpha = 0.5)
  g2 = ggplot(gd, aes(x=DESeq2, y=IDEAS_DCA)) + geom_point(alpha = 0.5)
  
  pdf(sprintf("figures/4_GSEA_p_val_compare_%s.pdf", ct1), 
      width=6, height=3)
  print(ggarrange(g1, g2, ncol = 2, nrow = 1))
  dev.off()
}

plot_gsea(gseaL, "CD8+Tcells_1")
#plot_gsea(gseaL, "PFC_Microglia")

# -------------------------------------------------------------------
# plot GSEA results for CD8+Tcells_1
# -------------------------------------------------------------------

ct1 = "CD8+Tcells_1"

#IDEAS_DCA = gesaL[[ct1]]$PS_dca_direct_Was
#summary(IDEAS_DCA$padj)
#IDEAS_DCA = IDEAS_DCA[which(IDEAS_DCA$padj < 0.05),]
#IDEAS_DCA$pathway = gsub("REACTOME_", "", IDEAS_DCA$pathway)

#df_gsea = as.data.frame(IDEAS_DCA)
#df_gsea = df_gsea[order(df_gsea$pval),]
#df_gsea$pathway = factor(df_gsea$pathway, levels=rev(df_gsea$pathway))
#p1 = ggplot(data=df_gsea, aes(x=pathway, y=-log10(pval))) +
#  geom_bar(stat="identity", fill="#56B4E9") + coord_flip() +
#  theme(legend.position = "none") + xlab("")

#ggsave("figures/2_GSEA_CD8+Tcells_1_IDEAS_DCA.pdf", p1, width=6, height=3)

for (i in 1:length(methods)){
  m = methods[i]
  cur_method = gesaL[[ct1]][[i]]
  cur_padj = cur_method$padj
  cat(paste(m, "  minimal adjusted pvalue: ",
            as.character(min(cur_padj, na.rm = TRUE)), "\n", sep = ""))
  
  row_sign_padj = which(cur_padj < 0.05)
  cat(paste("number of pathways with adjusted pvalue<0.05: ", 
            as.character(length(row_sign_padj)), "\n" ))

  if (length(row_sign_padj)>0){
    
    cur_method = cur_method[row_sign_padj,]
    cur_method$pathway = gsub("REACTOME_", "", cur_method$pathway)
    
    df_gsea = as.data.frame(cur_method)
    df_gsea = df_gsea[order(df_gsea$pval),]
    df_gsea$pathway = factor(df_gsea$pathway, levels=rev(df_gsea$pathway))
    p1 = ggplot(data=df_gsea, aes(x=pathway, y=-log10(pval))) +
      geom_bar(stat="identity", fill="#56B4E9") + coord_flip() +
      theme(legend.position = "none") + xlab("")
    
    cur_figurename = 
      paste("figures/4_GSEA_CD8+Tcells_1_", m, ".pdf", sep = "")
    ggsave(cur_figurename, p1, device = "pdf", width=12, 
           height=round(length(row_sign_padj)/4)+1)
  }
}




# -------------------------------------------------------------------
# scatter plot
# -------------------------------------------------------------------

pdf("figures/4_compare_pval_DESeq2_vs_IDEAS.pdf", width=4.5, height=4.5)

for(ct1 in cell_types){
  df1 = pv_list[[ct1]][,c("gene", "DESeq2", "PS_nb_Was", "PS_dca_direct_Was")]
  
  nms = c("DESeq2", "PS_nb_Was", "PS_dca_direct_Was")
  fun1 = function(x){-log10(x)}
  df2  = df1[, lapply(.SD, fun1), .SDcols = nms]
  df2[,gene:=df1$gene]
  setnames(df2, c("PS_nb_Was", "PS_dca_direct_Was"), 
           c("IDEAS_NB", "IDEAS_DCA"))
  
  df2[,label_NB:=""]
  df2[DESeq2 > 3 & IDEAS_NB < 1, label_NB:=gene]
  df2[DESeq2 < 1 & IDEAS_NB > 3, label_NB:=gene]
  
  df2[,color_NB:="g11"]
  df2[DESeq2 <= 2 & IDEAS_NB <= 2,  color_NB:="g00"]
  df2[DESeq2 > 3 & IDEAS_NB < 1, color_NB:="g01"]
  df2[DESeq2 < 1 & IDEAS_NB > 3, color_NB:="g01"]
  
  df2[,label_DCA:=""]
  df2[DESeq2 > 3 & IDEAS_DCA < 1, label_DCA:=gene]
  df2[DESeq2 < 1 & IDEAS_DCA > 3, label_DCA:=gene]
  
  df2[,color_DCA:="g11"]
  df2[DESeq2 <= 2 & IDEAS_DCA <= 2,  color_DCA:="g00"]
  df2[DESeq2 > 3 & IDEAS_DCA < 1, color_DCA:="g01"]
  df2[DESeq2 < 1 & IDEAS_DCA > 3, color_DCA:="g01"]
  
  df2$color_NB = factor(df2$color_NB, levels = c("g00", "g11", "g01", "g10"))
  df2$color_DCA = factor(df2$color_DCA, levels = c("g00", "g11", "g01", "g10"))
  
  g_NB = ggplot(df2, aes(x=DESeq2, y=IDEAS_NB, col=color_NB)) + 
    geom_point(size = 0.4) +
    scale_color_manual(values=c("#CCCCCC", "black", "#FC4E07", "#FC4E07")) +
    theme(legend.position = "none") + ggtitle(ct1)
  
  g_NB = g_NB + geom_text_repel(aes(label = df2$label_NB), size = 3.5, 
                                min.segment.length = 0, 
                                box.padding = unit(0.25, "lines"),
                                point.padding = unit(0.25, "lines")) 
  
  g_DCA = ggplot(df2, aes(x=DESeq2, y=IDEAS_DCA, col=color_DCA)) + 
    geom_point(size = 0.4) +  
    scale_color_manual(values=c("#CCCCCC", "black", "#FC4E07", "#FC4E07")) +
    theme(legend.position = "none") + ggtitle(ct1)
  
  g_DCA = g_DCA + geom_text_repel(aes(label = df2$label_DCA), size = 3.5, 
                                  min.segment.length = 0, 
                                  box.padding = unit(0.25, "lines"),
                                  point.padding = unit(0.25, "lines"))
  
  print(g_NB)
  print(g_DCA)
  
}

dev.off()






gc()

sessionInfo()
q(save="no")
