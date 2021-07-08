
#########################################################################
#                                                                       #
#                                                                       #
#             Part II: Defferential Expression Analysis                 #
#                                                                       #
#                                                                       #
#########################################################################
# once simulated genes, we will do the Differential Expression Analysis
# Here we will implement our method, plus the DESeq2 and the MAST analysis.

args = commandArgs(trailingOnly=TRUE)
args

if (length(args) < 5) {
  message("no enough arguments, using default values")
  r_mean   = 1.2     # The expected fold-changes in mean
  r_var    = 1.5     # The expected fold-changes in variances
  ncase    = 13      # case individuals
  nctrl    = 10      # control individuals
  ncell    = 360    # numbers of cells collected from each individuals.
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

if(ncell == 0){
  UNEQ_N_CELL = TRUE
}else{
  UNEQ_N_CELL = FALSE
}

if(UNEQ_N_CELL){
  config = sprintf("ncase_%d_nctrl_%d_unequal_n_cell", ncase, nctrl)
}else{
  config = sprintf("ncase_%d_nctrl_%d_ncell_%d", ncase, nctrl, ncell)
}

config = sprintf("%s_fold_mean_%.1f_var_%.1f", config, r_mean, r_var)
config

# ---------------------------------------------------------------
# initial setup
# ---------------------------------------------------------------

library(SAVER)

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

# ---------------------------------------------------------------
# load data
# ---------------------------------------------------------------

sim_data     = readRDS(sprintf("data/sim_data_%s.rds", config))
count_matrix = sim_data$count_matrix

# ---------------------------------------------------------------
# saver 
# ---------------------------------------------------------------

count_matrix_saver = saver(count_matrix,estimates.only = TRUE)

sim_data$count_matrix_saver = count_matrix_saver

# ---------------------------------------------------------------
# save results
# ---------------------------------------------------------------

saveRDS(dat_list, file=sprintf("data/sim_data_saver_%s.rds", config))

sessionInfo()

mem_used()
gc()

q(save = "no")
