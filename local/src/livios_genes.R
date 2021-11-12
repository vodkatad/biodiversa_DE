### from all genes by de rule obtain only livio's selected genes

library(tidyverse)

de_f <- snakemake@input[["tsv"]]
livio_genes <- snakemake@output[["de_genes"]]
de_excel <- snakemake@output[["excel"]]

#de_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici30/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"

de <- read.table(de_f, quote = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

de <- tibble::rownames_to_column(de, var = "genes")
de <- de %>% mutate(genes = gsub("H_", "", genes)) 


genes_livio <- c("RAD51", "RAD51C", "FBXO18", "SUMO1", "UBE2I", "BRCA1", "RRP1", 
                 "FBXO5", "RING1", "RFWD3", "UCHL3", "PARPBP", "BLM")  

de_fin <- de[de$genes %in% genes_livio,]

pvals <- de_fin$pvalue

de_fin$padjnew <- p.adjust(pvals, method="BH")
de_fin$padj <- NULL

names(de_fin)[names(de_fin) == 'padjnew'] <- 'padj'

write.table(de_fin, file = livio_genes, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

openxlsx::write.xlsx(de_fin, file = de_excel, sheetName = "Sheet1", 
                     col.names = TRUE, row.names = FALSE, append = FALSE)