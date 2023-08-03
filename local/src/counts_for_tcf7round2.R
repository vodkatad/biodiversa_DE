## samples data for rnaseq tcf7l2 2nd round

library(tidyverse)
library(readxl)

counts_f <- snakemake@input[["counts_or"]]
sample_f <- snakemake@input[["samples_or"]]
hmat_f_sh <- snakemake@output[["hmat_sh"]]
hmat_f_crc <- snakemake@output[["hmat_crc"]]

save.image("pippo.Rdata")
#counts_f <- "/mnt/cold1/bioinfotree/prj/ngs_standard_processing/dataset/TCF7L2_2nd/GEP.count.gz"
counts <- read.table(counts_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
counts <- as.data.frame(t(counts))

#samples <- read_xlsx("/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/TCF7L2_2nd/samples_names.xlsx")
samples <- read_xlsx(sample_f)
samples$...4 <- NULL
samples$...5 <- NULL
samples$sample <- gsub(" ", "_", samples$sample)

counts$sample <- rownames(counts)
counts <- merge(counts, samples, by="sample")
rownames(counts) <- counts$REPLICATES
counts$sample <- NULL
counts$`tube name` <- NULL
counts$REPLICATES <- NULL
counts$type <- substr(rownames(counts), 1,2)
counts_sh <- counts %>% filter(type == "sh")
counts_sh$type <- NULL
counts_crc <- counts %>% filter(!type== "sh")
counts_crc$type <- NULL

counts_sh <- as.data.frame(t(counts_sh))
counts_crc <- as.data.frame(t(counts_crc))
#res <- count

write.table(counts_sh, file=hmat_f_sh, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(counts_crc, file=hmat_f_crc, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
