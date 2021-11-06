# compared to 5f_gsea.R
# this step orders genes by the mean expression level


# compared to 3f, this code only deals with the results of 
# 5d_dca_direct.R

# this code is modified from 3f_gsea.R


# compared to 1f, the results delt with in this version consider more genes
# by setting filtering criterion to keep genes
# appearning in at last 10% of the cells




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


#pvals = fread(sprintf("res/3e_pvals_%s.tsv", grp))
pvals = fread(sprintf("res/5d_dca_direct_pvals_%s.tsv", grp))
pvals = pvals[, c(1, 4)]
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
# read in mean adjusted expression level 
# ------------------------------------------------------------------------

df_express = read.csv("res/5f_mean_expressions.csv", header = TRUE)

# check whether the genes match
# by checking the dca_direct pvalues
summary(abs(df_express$dca_direct + log10(pvals$PS_dca_direct_Was)))



# ------------------------------------------------------------------------
# GESA
# ------------------------------------------------------------------------


methods = names(df_express)[3:4]

gsea = list()

for(m1 in methods){
  #stats = unlist(-log10(pvals[[m1]]))
  stats = df_express[[m1]]
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

saveRDS(gsea, file = sprintf("res/5f_gsea_expression_%s.rds",grp))




gc()
sessionInfo()
q(save="no")
