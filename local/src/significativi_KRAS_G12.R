controllo <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC1598/M_C_vs_NT/CRC1598_treat_cutoff0.05-M_C.vs.NT.deseq2.tsv"
controllo <- read.table(controllo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
controllo <- controllo %>% filter(padj < 0.05)
controllo <- controllo %>% filter(abs(log2FoldChange) > 0.5849625)
controllo$genes <- rownames(controllo)

keys <- as.character(controllo$genes)
merge <- function(x) { paste(x, collapse=',') }
geneDescriptions <- mapIds(
  org.Hs.eg.db,
  keys = keys,
  column = "GENENAME",  
  keytype = "SYMBOL", 
  multiVals = merge
)
geneDescriptions[is.na(geneDescriptions) | geneDescriptions == "NA"] <- ""
controllo$description <- geneDescriptions
controllo <- controllo %>%
  mutate(description = geneDescriptions) %>%
  select(genes, description, everything())

controllo <- controllo[order(controllo$padj),]

write.xlsx(controllo, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC1598/M_C_vs_NT/CRC1598_treat_cutoff0.05-M_C.vs.NT.deseq2_sign_genedesc.xlsx",
           colNames = TRUE)

caso <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_vs_NT/CRC0031_treat_cutoff0.05-M_C.vs.NT.deseq2.tsv"
caso <- read.table(caso, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
caso <- caso %>% filter(padj < 0.05)
caso <- caso %>% filter(abs(log2FoldChange) > 0.5849625)
caso$genes <- rownames(caso)

keys <- as.character(caso$genes)
merge <- function(x) { paste(x, collapse=',') }
geneDescriptions <- mapIds(
  org.Hs.eg.db,
  keys = keys,
  column = "GENENAME",  
  keytype = "SYMBOL", 
  multiVals = merge
)
geneDescriptions[is.na(geneDescriptions) | geneDescriptions == "NA"] <- ""
caso$description <- geneDescriptions
caso <- caso %>%
  mutate(description = geneDescriptions) %>%
  select(genes, description, everything())
caso <- caso[order(caso$padj),]

write.xlsx(caso, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_vs_NT/CRC0031_treat_cutoff0.05-M_C.vs.NT.deseq2_sign_genedesc.xlsx",
           colNames = TRUE)

caso <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_B_vs_NT/CRC0031_treat_cutoff0.05-M_C_B.vs.NT.deseq2.tsv"
caso <- read.table(caso, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
caso <- caso %>% filter(padj < 0.05)
caso <- caso %>% filter(abs(log2FoldChange) > 0.5849625)
caso$genes <- rownames(caso)

keys <- as.character(caso$genes)
merge <- function(x) { paste(x, collapse=',') }
geneDescriptions <- mapIds(
  org.Hs.eg.db,
  keys = keys,
  column = "GENENAME",  
  keytype = "SYMBOL", 
  multiVals = merge
)
geneDescriptions[is.na(geneDescriptions) | geneDescriptions == "NA"] <- ""
caso$description <- geneDescriptions
caso <- caso %>%
  mutate(description = geneDescriptions) %>%
  select(genes, description, everything())
caso <- caso[order(caso$padj),]

write.xlsx(caso, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_B_vs_NT/CRC0031_treat_cutoff0.05-M_C_B.vs.NT.deseq2_sign_genedesc.xlsx",
           colNames = TRUE)

caso <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_T_vs_NT/CRC0031_treat_cutoff0.05-M_C_T.vs.NT.deseq2.tsv"
caso <- read.table(caso, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
caso <- caso %>% filter(padj < 0.05)
caso <- caso %>% filter(abs(log2FoldChange) > 0.5849625)
caso$genes <- rownames(caso)

keys <- as.character(caso$genes)
merge <- function(x) { paste(x, collapse=',') }
geneDescriptions <- mapIds(
  org.Hs.eg.db,
  keys = keys,
  column = "GENENAME",  
  keytype = "SYMBOL", 
  multiVals = merge
)
geneDescriptions[is.na(geneDescriptions) | geneDescriptions == "NA"] <- ""
caso$description <- geneDescriptions
caso <- caso %>%
  mutate(description = geneDescriptions) %>%
  select(genes, description, everything())
caso <- caso[order(caso$padj),]

write.xlsx(caso, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_T_vs_NT/CRC0031_treat_cutoff0.05-M_C_T.vs.NT.deseq2_sign_genedesc.xlsx",
           colNames = TRUE)

gsea_c <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C.vs.NT.tsv"
gsea_c <- read.table(gsea_c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_c$ID <- NULL
gsea_c <- gsea_c[order(gsea_c$p.adjust),]

gsea_c_t <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_T_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C_T.vs.NT.tsv"
gsea_c_t <- read.table(gsea_c_t, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_c_t$ID <- NULL
gsea_c_t <- gsea_c_t[order(gsea_c_t$p.adjust),]

gsea_c_a <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_B_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C_B.vs.NT.tsv"
gsea_c_a <- read.table(gsea_c_a, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_c_a$ID <- NULL
gsea_c_a <- gsea_c_a[order(gsea_c_a$p.adjust),]


wb <- createWorkbook()

addWorksheet(wb, "CRC0031_MRTX_Cet_GSEA_H")
addWorksheet(wb, "CRC0031_MRTX_Cet_Tram_GSEA_H")
addWorksheet(wb, "CRC0031_MRTX_Cet_Alp_GSEA_H")

writeData(wb, "CRC0031_MRTX_Cet_GSEA_H", gsea_c, startCol = 1, rowNames = FALSE)
writeData(wb, "CRC0031_MRTX_Cet_Tram_GSEA_H", gsea_c_t, startCol = 1, rowNames = FALSE)
writeData(wb, "CRC0031_MRTX_Cet_Alp_GSEA_H", gsea_c_a, startCol = 1, rowNames = FALSE)

saveWorkbook(wb, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/GSEA.xlsx", overwrite = TRUE)


