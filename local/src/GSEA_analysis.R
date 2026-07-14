### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
GSEA_r <- snakemake@output[["GSEA_r"]]
GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
type <- snakemake@wildcards[["msign"]]
rdata <- snakemake@output[["rdata"]] ### added by mv on 23/01/25 to save image with given wc
gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens", collection = type) %>% 
### Warning message:
### The `category` argument of `msigdbr()` is deprecated as of msigdbr 10.0.0.
### ℹ Please use the `collection` argument instead. 
### dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici
  dplyr::select(gs_name, gene_symbol) ### altrimenti chiede gli id numerici



em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)

write.table(em@result, file = GSEA_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

ridgeplot(em, showCategory = 20)
ggsave(GSEA_ridgeplot, width = 300, height = 107, useDingbats=FALSE, units = "mm")

#save.image(rdata) ### scommentato da mv 23/01/25 per deg 5vs4
save.image('GSEA.Rdata')
