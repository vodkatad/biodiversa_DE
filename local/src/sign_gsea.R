## 0.1 sui qvalues

library(tidyverse)
library(WriteXLS)

gseah <- snakemake@input[["h"]]
gseac <- snakemake@input[["c2"]]
results <- snakemake@output[["res"]]
excel <- snakemake@output[["xls"]]

#gseah <- "/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/GSEA_results_H_type_cutoff0.05-4Q.vs.1Q.tsv"
gseah <- read.table(gseah, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gseah <- gseah %>% filter(qvalues < 0.1)

#gseac <- "/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/GSEA_results_C2_type_cutoff0.05-4Q.vs.1Q.tsv"
gseac <- read.table(gseac, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gseac <- gseac %>% filter(qvalues < 0.1)

gseasign <- rbind(gseah, gseac)
gseasign$ID <- NULL
gseasign$Description <- NULL

write.table(gseasign, results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
WriteXLS(gseasign, excel, col.names = TRUE, row.names = TRUE)

