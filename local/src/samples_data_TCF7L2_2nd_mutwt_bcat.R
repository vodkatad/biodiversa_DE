library(tidyverse)

samples_f <- snakemake@input[["samples_or"]] 
meta <- snakemake@output[["meta"]]
type <- snakemake@wildcards[['tipo']]

d <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/b2_general/samples_data"
d <- read.table(d, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
d <- d %>% filter(!geno == "t2")

mut <- c("CRC0148", "CRC1331", "CRC0399", "CRC1278", "CRC0277", "CRC0327",
         "CRC1729", "CRC0152", "CRC0196", "CRC0059", "CRC0065", "CRC0464",
         "CRC0316", "CRC1239")
mut <- paste0("sh_", mut)
rownames(d) <- d$id

for (i in rownames(d)) {
  if (d[i, "model"] %in% mut) {
    d[i, "mut"] <- "MUT"
  } else {
    d[i, "mut"] <- "WT"
  }
}

d <- d %>% filter(mut == type)

write.table(d, file=meta,quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)