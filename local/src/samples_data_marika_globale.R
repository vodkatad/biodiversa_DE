## samples data marika for global deg
library(tidyverse)

c <- snakemake@input[["counts"]]
of <- snakemake@output[["o"]]

#c <- "/mnt/cold1/bioinfotree/prj/ngs_standard_processing/dataset/RNAseq_Marika_Trametinib_03_2026/GEP.count.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

delete <- c("Geneid", "CRC0322_APCmin", "CRC0322_APCplus", "CRC0327_APCmin", "CRC0327_APCplus", "CRC1502_NT_1", "CRC1502_TRAM_1", "CRC1477_TRAM_2")

samples <- as.data.frame(matrix(ncol = 1, nrow = length(colnames(c))))
samples$V1 <- colnames(c)
samples <- samples %>% filter(!V1 %in% delete)
samples$smodel <- substr(samples$V1, 1, 7)
samples$treat <- str_extract(samples$V1, "(?<=_)[^_]+(?=_)")
samples$replicate <- str_extract(samples$V1, "[^_]+$")
rownames(samples) <- samples$V1
samples$V1 <- NULL
samples$replicate <- paste0("R", samples$replicate)
samples$replicate <- as.factor(samples$replicate)
samples$treat <- as.factor(samples$treat)
samples$smodel <- as.factor(samples$smodel)

write.table(samples, file=of, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)