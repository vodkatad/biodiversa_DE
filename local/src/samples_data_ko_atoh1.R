library(stringr)

meta_f <- snakemake@input[["meta"]]
res <- snakemake@output[["tsv"]]

#ss <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/ko_atoh1/ss_ko_atoh1.txt"
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

write.table(ss, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
