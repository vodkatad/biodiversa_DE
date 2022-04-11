library(tidyverse)

input <- snakemake@input[["vsd"]]
output <- snakemake@output[["vsd_f"]]

#vsd <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(input, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
vsd <- vsd %>% remove_rownames %>% column_to_rownames(var="genes")

write.table(vsd, output, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
