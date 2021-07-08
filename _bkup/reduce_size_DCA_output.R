
library("data.table")

# -----------------------------------------------------------------
# try to reduce the output file of DCA using signif function
# -----------------------------------------------------------------

setwd("/fh/fast/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM3k10/")
t_mean       = fread("mean.tsv")
t_dispersion = fread("dispersion.tsv")
t_dropout    = fread("dropout.tsv")

rnms = as.character(unlist(t_mean[, 1]))
rnms[1:5]
length(rnms)
length(unique(rnms))

t_mean[1:2,1:5]
t_mean = data.matrix(t_mean[,-1])
t_mean = signif(t_mean, 4)
rownames(t_mean) = rnms
dim(t_mean)
t_mean[1:2,1:5]

write.table(t_mean, file = "mean_signif4.tsv", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

t_dispersion[1:2,1:5]
table(t_dispersion[,1] == rnms)
t_dispersion = data.matrix(t_dispersion[,-1])
t_dispersion = signif(t_dispersion, 4)
rownames(t_dispersion) = rnms
colnames(t_dispersion) = colnames(t_mean)
dim(t_dispersion)
t_dispersion[1:2,1:5]

write.table(t_dispersion, file = "dispersion_signif4.tsv", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

t_dropout[1:2,1:5]
table(t_dropout[,1] == rnms)
t_dropout = data.matrix(t_dropout[,-1])
t_dropout = signif(t_dropout, 4)
rownames(t_dropout) = rnms
colnames(t_dropout) = colnames(t_mean)
dim(t_dropout)
t_dropout[1:2,1:5]

write.table(t_dropout, file = "dropout_signif4.tsv", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

sessionInfo()
q(save="no")
