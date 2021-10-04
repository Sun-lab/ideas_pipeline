
# ========================================================================
# take arguments
# ========================================================================

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'CD8+Tcells_1' as default.\n")
  grp = "CD8+Tcells_1"
}else{
  eval(parse(text=args[[1]]))
}

grp


# ========================================================================
# libraries and path
# ========================================================================

library(MASS)
library(data.table)
library(doParallel)
library(doRNG)
library(qvalue)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(fgsea)
library(stringr)

# number of cores for multi-core computation
nCore = 4
# nCore = Sys.getenv("SLURM_CPUS_ON_NODE")
registerDoParallel(cores=nCore)
options(mc.cores=nCore)

RNGkind("L'Ecuyer-CMRG")
# ------------------------------------------------------------------------
# read in pathway information
# ------------------------------------------------------------------------

gmtfile_go_bp     = "../Autism/data/c5.bp.v7.1.symbols.gmt"
gmtfile_reactome  = "../Autism/data/c2.cp.reactome.v7.1.symbols.gmt"
pathways_go_bp    = gmtPathways(gmtfile_go_bp)
pathways_reactome = gmtPathways(gmtfile_reactome)

length(pathways_reactome)
summary(sapply(pathways_reactome, length))
reactome_genes = unique(unlist(pathways_reactome))
length(reactome_genes)

# ------------------------------------------------------------------------
# read in p-values
# ------------------------------------------------------------------------


pvals = fread(sprintf("res/1e_pvals_%s.tsv", grp))


dim(pvals)
head(pvals)

w_ens = grep("_", pvals$gene)
length(w_ens)

#pvals$gene[w_ens[1:20]]
genes = pvals$gene
#genes[w_ens] = str_extract(genes[w_ens], '(\\S+)(?=_)')  # gene symbol
#genes[w_ens][1:20]

pvals$gene_symbol = genes

# ------------------------------------------------------------------------
# filter genes in reactome annoations
# ------------------------------------------------------------------------

summary(sapply(pathways_reactome, length))
for(p1 in names(pathways_reactome)){
  pathways_reactome[[p1]] = intersect(pathways_reactome[[p1]], genes)
}
summary(sapply(pathways_reactome, length))


# ------------------------------------------------------------------------
# GESA
# ------------------------------------------------------------------------


methods = names(pvals)[2:17]

gsea = list()

for(m1 in methods){
  stats = unlist(-log10(pvals[[m1]]))
  names(stats) =  pvals$gene_symbol
  length(stats)
  
  stats = stats[unique(names(stats))]
  length(stats)
  stats = na.omit(stats)
  length(stats)
  
  set.seed(918)
  fgseaRes = fgseaMultilevel(pathways_reactome, stats, 
                             minSize=10, maxSize=1000)
  od1 = order(fgseaRes[,"padj"], -fgseaRes[,"NES"])
  fgseaRes   = fgseaRes[od1,]
  gsea[[m1]] = fgseaRes
}

lapply(gsea, head, n=10)
lapply(gsea, dim)

#lapply(gsea, function(x){x[which(x$pathway=="SFARI"),]})

saveRDS(gsea, file = sprintf("res/1f_gsea_%s.rds",grp))


# ------------------------------------------------------------------------
# compare GESA results across methods
# ------------------------------------------------------------------------

ps2 = merge(gsea$PS_nb_Was, gsea$PS_dca_direct_Was, by=c("pathway"), 
            suffixes = c( "_PS_nb_Was", "_PS_dca_direct_Was"))
ps2 = merge(gsea$DESeq2, ps2, by=c("pathway"))
change_labels = names(ps2)[2:8]
change_labels = paste(change_labels, "_DESeq2", sep = "")
names(ps2)[2:8] = change_labels

dim(ps2)
ps2[1:2,]

pvals.df = ps2[,.(pval_DESeq2, pval_PS_nb_Was, pval_PS_dca_direct_Was)]
dim(pvals.df)
pvals.df[1:2,]

colSums(is.na(pvals.df))

cor(-log10(pvals.df))
cor(-log10(pvals.df), method="spearman")

pvals.df.nona = na.omit(pvals.df)

cor(-log10(pvals.df.nona))
cor(-log10(pvals.df.nona), method="spearman")


gc()
sessionInfo()
q(save="no")
