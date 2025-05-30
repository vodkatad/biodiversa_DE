## deseq2 with chip samples with double crosslinking

library(tidyverse)

s <- snakemake@input[["samples"]]
tsv <- snakemake@output[["meta"]]

s <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/general/table_cases.tsv"
s <- read.table(s, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

chip <- c("CRC0148", "CRC0542", "CRC0291")

s <- s %>% filter(model %in% chip)
names(s)[names(s)=="geno"] <- "geno_chip"

write.table(s, tsv, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)