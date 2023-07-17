library(tidyverse)
library(readxl)

samples_f <- snakemake@input[["samples_or"]] 
meta <- snakemake@output[["meta"]]
type <- snakemake@wildcards[['geno']]

#samples_f <- read_xlsx("/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/TCF7L2_2nd/samples_names.xlsx")
samples <- read_xlsx(samples_f)
samples$...4 <- NULL
samples$...5 <- NULL
samples <- samples[-c(133:171),]
samples$id <- samples$REPLICATES
samples$model <- substr(samples$REPLICATES, 1,7)
samples$geno <- substr(samples$REPLICATES, 9, 10)
samples$replicates <- substr(samples$REPLICATES, 12, 13)
samples$sample <- NULL
samples$`tube name` <- NULL
samples$REPLICATES <- NULL

samples <- samples %>% filter(geno == type)

mut <- c("CRC0148", "CRC1331", "CRC0399", "CRC1278", "CRC0277", "CRC0327",
         "CRC1729", "CRC0152", "CRC0196", "CRC0059", "CRC0065", "CRC0464",
         "CRC0316", "CRC1239")

samples <- as.data.frame(samples)
rownames(samples) <- samples$id

for (i in rownames(samples)) {
  if (samples[i, "model"] %in% mut) {
    samples[i, "mut"] <- "MUT"
  } else {
    samples[i, "mut"] <- "WT"
  }
}

samples$geno_mut <- paste0(samples$geno, "_", samples$mut)
rownames(samples) <- NULL
samples$geno <- NULL
samples$mut <- NULL

#samples <- samples %>% filter(model == case)
write.table(samples, file=meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)