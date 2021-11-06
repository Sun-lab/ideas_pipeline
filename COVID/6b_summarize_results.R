# compared with 6a, this version only deals with the results
# based on 5f_gsea_expression.R



# compared with 4a, this version only deals with the results
# based on 4d_dca_direct PS_dca_direct_Was


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

methods = c("paper", "email")


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

sink("res/6b_GSEA.txt")

for(ct1 in cell_types){
  
  cat("\n---------------------------------------------------\n")
  cat(ct1)
  cat("\n---------------------------------------------------\n")
  
  gsea = readRDS(sprintf("res/5f_gsea_expression_%s.rds", ct1))
  gsea = gsea[methods]
  
  x1   = lapply(gsea, f1)
  x1[sapply(x1, is.null)] = NULL
  
  if(length(x1) > 0){ print(x1) }
  
  gesaL[[ct1]] = gsea
}
sink()



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
      paste("figures/6b_GSEA_CD8+Tcells_1_", m, "_expression.pdf", sep = "")
    ggsave(cur_figurename, p1, device = "pdf", width=12, 
           height=round(length(row_sign_padj)/4)+1, limitsize = FALSE)
    
    df_gsea$leadingEdge <- 
      vapply(df_gsea$leadingEdge, paste, collapse = ", ", character(1L))
    
    write.csv(df_gsea, 
              file = sprintf("res/6b_GSEA_sign_pathways_CD8+Tcells_1_%s.csv", m), 
              row.names = FALSE)
  }
}








gc()

sessionInfo()
q(save="no")
