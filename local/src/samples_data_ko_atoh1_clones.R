library(stringr)
library(tidyverse)

meta_f <- snakemake@input[["meta"]]
res <- snakemake@output[["tsv"]]

#meta_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/ko_atoh1/ss_ko_atoh1.txt"
ss <- read.table(meta_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ss$V1 <- NULL                 
ss$model <- substr(ss$V2, 1, 7)
ss$V3 <- gsub("NO_EGF", "NO.EGF", ss$V2)
ss$V3 <- gsub("CAS9_1", "CAS9.1", ss$V3)
ss$V3 <- gsub("CAS9_2", "CAS9.2", ss$V3)


ss[c("model", "geno", "trattamento")] <- str_split_fixed(ss$V3, "_",3)
ss$V2 <- NULL
rownames(ss) <- ss$V3
ss$V3 <- NULL

rownames(ss) <- gsub("NO.EGF", "NO_EGF", rownames(ss))
rownames(ss) <- gsub("CAS9.1", "CAS9_1", rownames(ss))
rownames(ss) <- gsub("CAS9.2", "CAS9_2", rownames(ss))

ss$trattamento <- gsub("NO.EGF", "NO_EGF", ss$trattamento)
ss$geno <- gsub("CAS9.1", "CAS9_1", ss$geno)
ss$geno <- gsub("CAS9.2", "CAS9_2", ss$geno)

#caso <- c("CRC0327_CAS9_2_NO_EGF", "CRC0327_CAS9_2_CTX", "CRC0327_CAS9_1_NO_EGF", "CRC0327_CAS9_1_CTX")
caso <- c("CRC0322_CL7_NO_EGF", "CRC0322_CL7_CTX", "CRC0322_CL8_NO_EGF", "CRC0322_CL8_CTX",'CRC0322_CL11_NO_EGF','CRC0322_CL11_CTX')

ss$V3 <- rownames(ss)
ss <- ss %>% filter(V3 %in% caso)
ss$V3 <- NULL
# ss$replicati <- substr(ss$geno, 5,6)
# ss$replicati <- gsub("_", "R", ss$replicati)
# ss$geno <- gsub("_1", "", ss$geno)
# ss$geno <- gsub("_2", "", ss$geno)

write.table(ss, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)