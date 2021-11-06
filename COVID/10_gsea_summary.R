# check the overlap of singifcant pathways 

library(readxl)

fnm = "res/COVID_GSEA_sign_pathways_CD8+Tcells_1.xlsx"
excel_sheets(fnm)

DESeq2 = read_excel(fnm, sheet="DESeq2")
dim(DESeq2)

ideas = read_excel(fnm, sheet="PS_nb_Was")
dim(ideas)

dca = read_excel(fnm, sheet="PS_dca_direct_Was")
dim(dca)

saver = read_excel(fnm, sheet="PS_saver_direct_Was")
dim(saver)

meanE = read_excel(fnm, sheet="mean_expression")
dim(meanE)

meanE[1:2,]
summary(meanE$padj)

table(DESeq2$pathway %in% meanE$pathway)
table(ideas$pathway %in% meanE$pathway)
table(dca$pathway   %in% meanE$pathway)
table(saver$pathway %in% meanE$pathway)

options(width = 300)


ideas_E = merge(ideas[,c("pathway", "padj", "pval")], meanE[,c("pathway", "padj")], 
                by=c("pathway"), all.x = TRUE)
ideas_E[order(ideas_E$pval),]


dca_E = merge(dca[,c("pathway", "padj", "pval")], meanE[,c("pathway", "padj")], 
                by=c("pathway"), all.x = TRUE)
dca_E[order(dca_E$pval),]

saver_E = merge(saver[,c("pathway", "padj", "pval")], meanE[,c("pathway", "padj")], 
              by=c("pathway"), all.x = TRUE)
saver_E[order(saver_E$pval),]

gc()
sessionInfo()
q(save="no")
