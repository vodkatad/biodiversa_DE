library(stringr)
library(tidyverse)

meta_f <- snakemake@input[["meta"]]
res <- snakemake@output[["tsv"]]
#meta_f<-'/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/bulkRNA_pairedscRNA_KRASEGF/metadata_bulkRNA_pairedscRNA_KRASEGF.txt'


ss <- read.table(meta_f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ss$V1 <- NULL                 
ss$model <- substr(ss$V2, 1, 7)
ss$time <- gsub(".*_(\\d+)d_.*", "\\1", ss$V2)

rownames(ss)<-ss$V2
ss$V2<-NULL

write.table(ss, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
