### GO enrichment analysis
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)


gene_res_f <- snakemake@input[["geni"]]
gene_univ_f <- snakemake@input[["universo"]] 
plot_output <- snakemake@output[["plot_out"]]
print(plot_output)
tab_out <- snakemake@output[["tab_out"]]
print(tab_out)
universe<-read.table(gene_univ_f, quote = "", sep = "\t", header = TRUE)
fpkm<-as.data.frame(universe)
universe$gene<-row.names(universe)


gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE,row.names = 1)


geneList <- as.character(gene_res_df$gene)


geneUni <- universe$gene
geneUni <- as.character(geneUni)


egocc <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "CC",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

egomf <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "MF",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

egobp <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)


pdf(plot_output)
barplot(egocc, showCategory = 20,fontsize_row = 5,fontsize_col = 5)+ggtitle("CC")

barplot(egomf, showCategory = 20,fontsize_row = 5,fontsize_col = 5)+ggtitle("MF")

barplot(egobp, showCategory = 20,fontsize_row = 5,fontsize_col = 5)+ggtitle("BP")

graphics.off()

egocc@result$ontology <- "CC"
egobp@result$ontology <- "BP"
egomf@result$ontology <- "MF"
egoall_df <- rbind(egocc@result, egobp@result, egomf@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method='BH')
write.table(egoall_df, file = tab_out, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
