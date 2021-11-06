# this file combines the qvalue of DESeq2, dca_direct
# and the pvalues of log(mean), log(theta) regression 
# into one table

library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(transport)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)

library(grid)
library(gridExtra)
#big_font <- theme_grey(base_size =  16)

theme_set(theme_classic())


data.dir  = "./data"
data.dca.dir = "../../ideas_data/Autism/dca_PFC_all"

grp = "PFC_L2_3"
grp1 = "L2_3"



# -------------------------------------------------------------------
# load means and variances pvalues from log transformation
# and linear regression on covariates
# from formula and pmf based approach
# -------------------------------------------------------------------

dca_pvalues = read.csv(
  "res/step3b_formula_covariates_pvals_all.csv", 
  header = TRUE)
dim(dca_pvalues)
dca_pvalues[1:2, ]

theta_pvalues = read.csv(
  "res/step10b_formula_covariates_pvals_all.csv", 
  header = TRUE)
dim(theta_pvalues)
theta_pvalues[1:2, ]

dca_pvalues$theta_formula = theta_pvalues$theta_formula
dim(dca_pvalues)
dca_pvalues[1:2, ]

# -------------------------------------------------------------------
# load DESeq2 and dca_direct qvalues 
# -------------------------------------------------------------------

cell_types = scan("cell_types.txt", what=character())
cell_types = sort(cell_types)
cell_types

i = 8
ct1 = cell_types[i]

methods = c("DESeq2", "rank_sum", "MAST_glm", "MAST_glmer", 
            "PS_nb_Was", "PS_dca_direct_Was", "PS_saver_direct_Was")

q1  =  fread(sprintf("res/step1l_qvals_%s.tsv", ct1))

table(dca_pvalues$gene == q1$gene)
table(theta_pvalues$gene == q1$gene)


df = cbind(q1$gene, q1$DESeq2, q1$PS_dca_direct_Was, 
           dca_pvalues$mean_formula, theta_pvalues$theta_formula)

colnames(df) = c("gene", "DESeq2_qvalue", "PS_dca_direct_Was_qvalue", 
                 "log_mean_pvalue", "log_pseudo_theta_pvalue")

write.csv(df, 
          file = "res/Autism_PFC_L2_3_DESeq2_dca_direct_q_log_mean_theta_p.csv", 
          row.names = FALSE)

