
library(ggcorrplot)
library(data.table)
library(ggpubr)
library(ggrepel)
library(reshape2)

theme_set(theme_classic())

options(width = 120)

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

pi0     = matrix(NA, nrow=length(cell_types), ncol=9)
ngenes  = rep(NA, length(cell_types))
or_list = fp_list = pv_list = list()

pdf("figures/step2_pval_hist.pdf", width=9, height=9)
par(mfrow=c(3,3), bty="n", mar=c(5,4,3,1))

for(i in 1:length(cell_types)){
  ct1 = cell_types[i]
  p1  = fread(sprintf("res/step1e_pvals_%s.tsv", ct1))
  pv_list[[ct1]]= p1
  ngenes[i] = nrow(p1)
  
  odds_ratio = fisher_pvl = matrix(NA, nrow=9, ncol=9)
  colnames(odds_ratio) = rownames(odds_ratio) = names(p1)[-1]
  colnames(fisher_pvl) = rownames(fisher_pvl) = names(p1)[-1]
  
  for(j in 2:10){
    
    hist(p1[[j]], main=paste0(names(p1)[j], "\n", ct1), 
         xlab="p-value", breaks = 20)
    pi0[i,j-1] = min(1, 2*mean(p1[[j]] > 0.5, na.rm=TRUE))
    
    pj = p1[[j]] < 0.05
    if(sum(pj, na.rm=TRUE) == 0){ next }
    
    for(k in 2:10){
      pk = p1[[k]] < 0.05
      if(sum(pk, na.rm=TRUE) == 0){ next }
      
      fjk = fisher.test(pj, pk)
      odds_ratio[j-1, k-1] = fjk$estimate
      fisher_pvl[j-1, k-1] = fjk$p.value
    }
  }
  
  or_list[[ct1]] = odds_ratio
  fp_list[[ct1]] = fisher_pvl
  
}

dev.off()


df1 = data.frame(cell_type = gsub("PFC_", "", cell_types), n_genes = ngenes)
df1$cell_type = factor(df1$cell_type, levels = df1$cell_type)

g2  = ggplot(data=df1, aes(x=cell_type, y=n_genes)) +
  geom_bar(stat="identity", fill="gray") + 
  ylab("number of genes") + xlab("") + 
  coord_flip() + theme_classic() 

ggsave("figures/step2_barplot_n_genes.pdf", g2, width=3, height=3.5)

pdf("figures/step2_odds_ratio.pdf", width=6, height=6)
par(mfrow=c(1,1), bty="n", mar=c(5,4,3,1))

for(ct1 in cell_types){
  or1 = or_list[[ct1]]
  or1[which(or1 > 10)] = 10
  g1 = ggcorrplot(or1, tl.srt=90) + ggtitle(gsub("PFC_", "", ct1)) +
    scale_fill_gradient2(limit = c(0,10), low = "blue", 
                         high =  "red", mid = "white", 
                         midpoint = 5)
  print(g1)
}
dev.off()

# -------------------------------------------------------------------
# summarize the proportion of non-nulls
# -------------------------------------------------------------------

colnames(pi0) = names(p1)[-1]
rownames(pi0) = cell_types

pi1 = t(1 - pi0)

pi1
rownames(pi1) = gsub("dca_direct", "dca", rownames(pi1))
summary(c(pi1))

g1 = ggcorrplot(pi1, tl.srt=90)
g1 = g1 + scale_fill_gradient2(limit = c(0,0.62), low = "blue", 
                               high =  "red", mid = "white", 
                               midpoint = 0.31)

ggsave("figures/step2_pi1_complete.pdf", g1, width=5, height=5)

pi1_sub = pi1[c("DESeq2", "PS_nb_Was", "PS_dca_Was"),]
rownames(pi1_sub) = c("DESeq2", "IDEAS_NB", "IDEAS_DCA")

g1 = ggcorrplot(pi1_sub, tl.srt=90)
g1 = g1 + scale_fill_gradient2(limit = c(0,0.62), low = "blue", 
                               high =  "red", mid = "white", 
                               midpoint = 0.31)

ggsave("figures/step2_pi1_PS_Was.pdf", g1, width=4, height=5)

pi1_sub = as.data.frame(round(t(pi1_sub),3))
rownames(pi1_sub) = gsub("PFC_", "", rownames(pi1_sub))

table(df1$cell_type == rownames(pi1_sub))
pi1_sub
df1 = cbind(df1, pi1_sub)
df1

write.table(df1, file="res/Fig4A_n_genes.txt", sep="\t", quote=FALSE, 
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

f2 <- function(x){
  x[which(x$pathway=="SFARI"), c("pval")]
}

f3 <- function(x){
  x[which(x$pathway=="SFARI"), c("NES")]
}

sfari_pval = sfari_NES = NULL

gesaL = list()

sink("res/step2_GSEA.txt")

for(ct1 in cell_types){
  
  cat("\n---------------------------------------------------\n")
  cat(ct1)
  cat("\n---------------------------------------------------\n")
  
  gsea = readRDS(sprintf("res/step1f_gsea_%s.rds", ct1))
  x1   = lapply(gsea, f1)
  x1[sapply(x1, is.null)] = NULL
  
  if(length(x1) > 0){ print(x1) }
  
  x2  = unlist(sapply(gsea, f2))
  sfari_pval = rbind(sfari_pval, x2)

  x3   = unlist(sapply(gsea, f3))
  sfari_NES = rbind(sfari_NES, x3)
  
  gesaL[[ct1]] = gsea
}
sink()

rownames(sfari_pval) = rownames(sfari_NES) = cell_types

signif(sfari_pval,2)
signif(sfari_NES,2)

cols = c("DESeq2.pval", "PS_nb_Was.pval", "PS_dca_direct_Was.pval")
w_ct = which(rowSums(sfari_pval[,cols] < 0.1) > 0)
signif(sfari_pval[w_ct, cols],2)
signif(sfari_NES[w_ct, match(cols, colnames(sfari_pval))],2)

# ------------------------------------------------------------------------
# compare p-values of SFARI genes
# ------------------------------------------------------------------------

fnm = "data/SFARI-Gene_genes_10-31-2019release_04-08-2020export.csv"
SFARI_table0 = read.csv(fnm, as.is=TRUE)
dim(SFARI_table0)
SFARI_table0[1:2,]

table(SFARI_table0$gene.score, useNA="ifany")
table(SFARI_table0$syndromic, useNA="ifany")
table(SFARI_table0$gene.score <= 3, SFARI_table0$syndromic, useNA="ifany")

sum(SFARI_table0$gene.score <= 3 | SFARI_table0$syndromic==1, na.rm = TRUE)

# the cell types that are associated with SFARI among at least 
# of the the three methods at p-value 0.1
cts_SFARI = rownames(sfari_pval)[w_ct]

df_all  = NULL
pi1_all = NULL
fun.pi0 = function(v){ 2 * length(which(v>0.5)) / length(v) }

for(ct1 in cts_SFARI){
  # ct1 = "PFC_L2_3"
  cat("\n-----------------------------------------------------\n")
  cat(ct1)
  cat("\n-----------------------------------------------------\n")
  
  pvals = pv_list[[ct1]]
  
  SFARI_table = SFARI_table0[SFARI_table0$gene.symbol %in% pvals$gene, ]
  dim(SFARI_table)
  
  table(SFARI_table$syndromic, SFARI_table$gene.score)
  w2kp = SFARI_table$syndromic == 1 | SFARI_table$gene.score %in% 1:3
  genes_in_SFARI = SFARI_table$gene.symbol[w2kp]
  length(genes_in_SFARI)
  genes_in_SFARI[1:5]
  
  wSFARI = which(pvals$gene %in% genes_in_SFARI)
  length(wSFARI)
  pvals$SFARI = rep("No", nrow(pvals))
  pvals$SFARI[wSFARI] = "Yes"
  
  pcut = 0.05
  f1 = fisher.test(pvals$SFARI, pvals$DESeq2 < pcut)
  f2 = fisher.test(pvals$SFARI, pvals$PS_nb_Was < pcut)
  f3 = fisher.test(pvals$SFARI, pvals$PS_dca_direct_Was  < pcut)
  
  print(f1)
  print(f2)
  print(f3)
  
  methods = c("DESeq2", "IDEAS_NB", "IDEAS_DCA")
  methods_idx = 3:1
  
  df = data.frame(cell_type=rep(ct1,3), methods_idx = methods_idx, methods=methods,
                  odds   = c(f1$estimate, f2$estimate, f3$estimate),
                  CILow  = c(f1$conf.int[1], f2$conf.int[1], f3$conf.int[1]), 
                  CIHigh = c(f1$conf.int[2], f2$conf.int[2], f3$conf.int[2]),
                  pval   = c(f1$p.value, f2$p.value, f3$p.value))

  p = ggplot(df, aes(x = odds, y = methods_idx)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, 
                   height = .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    theme(panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = methods_idx, labels = methods) +
    scale_x_continuous(breaks = seq(0,2,0.5) ) +
    coord_trans(x = "log10") + ylab("") + 
    xlab("Odds ratio for SFARI genes") + ggtitle(gsub("PFC_", "", ct1))
  
  ggsave(sprintf("figures/step2_SFARI_%s_odds_ratio.pdf", ct1), p, 
         width=3, height=1.3)
  
  df_all = rbind(df_all, df)
  
  k = 0
  methods_label = c("DESeq2", "PS_nb_Was", "PS_dca_direct_Was")
  pi0 = matrix(NA, nrow=length(methods), ncol=2)
  
  for(m1 in methods_label){
    k = k + 1
    pi0[k,] = tapply(pvals[[m1]], pvals$SFARI, fun.pi0)
  }
  
  colnames(pi0) = c("not_SFARI", "SFARI")
  pi1_df = data.frame(cell_type=rep(ct1,3), methods=methods, 1-pi0)
  pi1_all = rbind(pi1_all, pi1_df)
  pi1 = melt(pi1_df, id=c("cell_type","methods"), variable.name="group", 
             value.name="non_null_prop")
  
  pi1$methods = factor(pi1$methods, 
                       levels=c("DESeq2", "IDEAS_NB", "IDEAS_DCA"))
  p2 = ggplot(pi1, aes(x = methods, y = non_null_prop, fill=group)) + 
    geom_bar(stat="identity", position=position_dodge()) + 
    scale_fill_brewer(palette="Dark2") + ggtitle(gsub("PFC_", "", ct1)) + 
    theme(axis.title.x = element_blank()) + ylab("non-null proportion") + 
    theme(legend.position = "topleft")
    
  ggsave(sprintf("figures/step2_SFARI_%s_non_null.pdf", ct1), p2, 
         width=3, height=2.5)
  
}

df_all
pi1_all

df_all$cell_type  = gsub("PFC_", "", df_all$cell_type)
pi1_all$cell_type = gsub("PFC_", "", pi1_all$cell_type)

df_all$methods = factor(df_all$methods, 
                        levels=c("DESeq2", "IDEAS_NB", "IDEAS_DCA"))

pi1_all$methods = factor(pi1_all$methods, 
                        levels=c("DESeq2", "IDEAS_NB", "IDEAS_DCA"))

gb1 = ggplot(data=df_all, aes(x=cell_type, y=-log10(pval), fill=methods)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.25) + 
  annotate(geom = "text", y =1.4, x = 0.5, label ="pval=0.05", size = 4, hjust=0)

ggsave("figures/step2_SFARI_fisher_test.pdf", gb1, width=5, height=2)

gb2 = ggplot(data=pi1_all, aes(x=not_SFARI, y=SFARI, 
                               color=methods, shape=cell_type)) +
  geom_point() + geom_abline(intercept = 0, slope=1, linetype="dashed", size=0.25)

ggsave("figures/step2_SFARI_non_null_prop.pdf", gb2, width=6.5, height=4.5)

# -------------------------------------------------------------------
# check GSEA results for one cell type
# -------------------------------------------------------------------

plot_gsea <- function(gseaL, ct1){
  # ct1 = "PFC_Microglia"
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
  
  pdf(sprintf("figures/step2_GSEA_p_val_compare_%s.pdf", ct1), 
      width=6, height=3)
  print(ggarrange(g1, g2, ncol = 2, nrow = 1))
  dev.off()
}

plot_gsea(gseaL, "PFC_IN-SST")
plot_gsea(gseaL, "PFC_Microglia")

# -------------------------------------------------------------------
# plot GSEA results for PFC_Microglia
# -------------------------------------------------------------------

ct1 = "PFC_Microglia"
IDEAS_DCA = gesaL[[ct1]]$PS_dca_direct_Was
IDEAS_DCA = IDEAS_DCA[which(IDEAS_DCA$padj < 0.05),]
IDEAS_DCA$pathway = gsub("REACTOME_", "", IDEAS_DCA$pathway)

df_gsea = as.data.frame(IDEAS_DCA)
df_gsea = df_gsea[order(df_gsea$pval),]
df_gsea$pathway = factor(df_gsea$pathway, levels=rev(df_gsea$pathway))
p1 = ggplot(data=df_gsea, aes(x=pathway, y=-log10(pval))) +
  geom_bar(stat="identity", fill="#56B4E9") + coord_flip() +
  theme(legend.position = "none") + xlab("")

ggsave("figures/step2_GSEA_microglia_IDEAS_DCA.pdf", p1, width=6, height=3)

# -------------------------------------------------------------------
# scatter plot
# -------------------------------------------------------------------

pdf("figures/step2_compare_pval_DESeq2_vs_IDEAS.pdf", width=4.5, height=4.5)

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
    theme(legend.position = "none") + ggtitle(gsub("PFC_", "", ct1))
  
  g_NB = g_NB + geom_text_repel(aes(label = df2$label_NB), size = 3.5, 
                                min.segment.length = 0, 
                                box.padding = unit(0.25, "lines"),
                                point.padding = unit(0.25, "lines")) 
  
  g_DCA = ggplot(df2, aes(x=DESeq2, y=IDEAS_DCA, col=color_DCA)) + 
    geom_point(size = 0.4) +  
    scale_color_manual(values=c("#CCCCCC", "black", "#FC4E07", "#FC4E07")) +
    theme(legend.position = "none") + ggtitle(gsub("PFC_", "", ct1))
  
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
