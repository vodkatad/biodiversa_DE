## gsea results tutti i risultati

library(tidyverse)
library(openxlsx)

h <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_H_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
h <- read.table(h, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

c2 <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C2_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
c2 <- read.table(c2, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

c6 <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C6_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
c6 <- read.table(c6, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "GSEA_results_Hallmark")
addWorksheet(OUT, "GSEA_results_C2")
addWorksheet(OUT, "GSEA_results_C6")

# Write the data to the sheets
writeData(OUT, sheet = "GSEA_results_Hallmark", x = h)
writeData(OUT, sheet = "GSEA_results_C2", x = c2)
writeData(OUT, sheet = "GSEA_results_C6", x = c6)

# Export the file
saveWorkbook(OUT, "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_all.xlsx")
