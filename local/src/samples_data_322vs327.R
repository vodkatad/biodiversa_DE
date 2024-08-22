library(stringr)
library(tidyverse)

meta_f <- snakemake@input[["meta"]]
res <- snakemake@output[["tsv"]]

#meta_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/322_vs_327/ss_322_vs_327.txt"
ss <- read.table(meta_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ss$model <- substr(ss$V1, 5, 7)
ss$V2 <- gsub("NO_EGF", "NO.EGF", ss$V1)
ss$V2 <- gsub("CAS9_1", "CAS9.1", ss$V2)
ss$V2 <- gsub("CAS9_2", "CAS9.2", ss$V2)
ss$V2<-gsub("EGF_1","EGF.1",ss$V2)
ss$V2<-gsub("EGF_2","EGF.2",ss$V2)
ss$V2<-gsub("NOEGF_1","NO.EGF.1",ss$V2)
ss$V2<-gsub("NOEGF_2","NO.EGF.2",ss$V2)



ss[c("model", "geno", "trattamento")] <- str_split_fixed(ss$V2, "_",3)
ss$model <- substr(ss$V1, 5, 7)
ss$V2 <- NULL
rownames(ss) <- ss$V1


ss$trattamento <- gsub("EGF.1", "EGF", ss$trattamento)
ss$trattamento <- gsub("EGF.2", "EGF", ss$trattamento)
ss$trattamento <- gsub("NO.EGF.1", "NO_EGF", ss$trattamento)
ss$trattamento <- gsub("NO.EGF.2", "NO_EGF", ss$trattamento)
ss$trattamento<-gsub('NO.EGF','NO_EGF',ss$trattamento)
ss$geno <- gsub("CAS9.1", "CAS9", ss$geno)
ss$geno <- gsub("CAS9.2", "CAS9", ss$geno)
ss$batch<-ifelse(grepl("CRC", ss$V1), "new", "old")
ss$model<-paste0('CRC0',ss$model)

ss$V1<-NULL

write.table(ss, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)