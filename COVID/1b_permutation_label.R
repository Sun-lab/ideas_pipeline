#library(Matrix)
library(data.table)
#library(dplyr)

data.dir = "data"
# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"))
dim(cell_info)
cell_info[1:2,]

table(cell_info$region)

cell_info = cell_info[which(cell_info$region=="PFC"),]

#individual level cell info
cell_info_ind = unique(paste(cell_info$diagnosis, cell_info$sample, sep=":"))

#map from ind level info to cell level info
match_ind_index = match(cell_info$individual,unique(cell_info$individual))

case_index = grep("ASD",cell_info_ind)
length(case_index)

set.seed(1)
perm_table=matrix(ncol=100,nrow=length(match_ind_index))
dim(perm_table)

# try to get balanced case/control sample here, 
# so only choose 20 out of 23 samples for permutations.
for(i in 1:100){
  select_index = seq(1:23)[-case_index[sample.int(13,3)]]
  case_index_perm=sample.int(20,10)
  perm_label=rep(NA,23)
  perm_label[select_index]="Control"
  perm_label[select_index[case_index_perm]]="ASD"
  perm_table[,i]=perm_label[match_ind_index]
}

colnames(perm_table)=paste0("p",1:100)
##to use it ####### add current lines after 
#cell_info = cell_info[which(cell_info$region=="PFC"),]

#eg:
#cell_info$sample = perm_table[,i] #for ith permutation
#dat = dat[,!is.na(cell_info$sample)]
#cell_info = cell_info[!is.na(cell_info$sample),]

# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

fwrite(data.frame(perm_table), file=sprintf("data/step1b_permlabel.tsv"), 
       sep = "\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

gc()

sessionInfo()
q(save="no")
