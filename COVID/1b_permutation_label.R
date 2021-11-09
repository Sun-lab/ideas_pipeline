library(Matrix)
library(data.table)


data.dir = "../../ideas_data/COVID/PBMC_10x"

args=(commandArgs(TRUE))
args

if (length(args) != 1) {
  message("one argument is expected, use 'CD8+Tcells_1' as default.\n")
  grp = "CD8+Tcells_1"
}else{
  eval(parse(text=args[[1]]))
}

grp

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = fread(file.path(data.dir, "meta.tsv"))
dim(cell_info)
cell_info[1:2,]

table(cell_info$id.celltype)


# ------------------------------------------------------------------------
# read in count data of cell type
# ------------------------------------------------------------------------

dat = readRDS(file.path(data.dir, sprintf("ct_mtx/%s.rds", grp)))
dim(dat)
class(dat)
dat[1:5,1:4]


# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat) %in% cell_info$cell)
meta = cell_info[match(colnames(dat), cell_info$cell),]
dim(meta)
meta[1:2,]


# ------------------------------------------------------------------------
# filter out cells from control samples
# ------------------------------------------------------------------------

table(meta$donor, meta$group_per_sample)

meta_covid = meta[which(meta$group_per_sample != "control"),]
dim(meta_covid)

table(meta_covid$donor, meta_covid$group_per_sample)

df_donor = as.data.frame(table(meta_covid$donor))
donor2kp = df_donor$Var1[which(df_donor$Freq >= 10)]

meta2kp = meta_covid[which(meta_covid$donor %in% donor2kp),]
dim(meta2kp)
table(meta2kp$donor)
length(unique(meta2kp$donor))


#concatenate disease status and donor id
#cell_info_ind = unique(paste(meta2kp$group_per_sample, meta2kp$donor, sep=":"))
status_and_donor = unique(paste(meta2kp$group_per_sample, meta2kp$donor, sep=":"))
status_and_donor

# double check that the order in status_and_donor matches that
# from meta2kp$donor
extract_donor <- function(item){
  return(unlist(strsplit(item, "[:]"))[2])
}

donor2check = unlist(lapply(status_and_donor, extract_donor))
donor2check

table(donor2check == unique(meta2kp$donor))


#map from donor level info to cell level info
match_ind_index = match(meta2kp$donor,unique(meta2kp$donor))

case_index = grep("severe",status_and_donor)
length(case_index)

set.seed(1)
perm_table=matrix(ncol=10,nrow=length(match_ind_index))
dim(perm_table)

# try to get balanced case/control sample here, 
# so only choose 7 out of 10 samples for severe.
for(i in 1:10){
  select_index = seq(1:17)[-case_index[sample.int(10,3)]]
  case_index_perm=sample.int(14,7)
  perm_label=rep(NA,17)
  perm_label[select_index]="mild"
  perm_label[select_index[case_index_perm]]="severe"
  perm_table[,i]=perm_label[match_ind_index]
}

colnames(perm_table)=paste0("p",1:10)

# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

fwrite(data.frame(perm_table), 
       file=sprintf("../../ideas_data/COVID/1b_permlabel.tsv"), 
       sep = "\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

gc()

sessionInfo()
q(save="no")
