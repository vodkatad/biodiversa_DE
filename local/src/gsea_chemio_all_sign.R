### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
imagine <- snakemake@output[["rdata"]]
results <- snakemake@output[["tsv"]]
#GSEA_r <- snakemake@output[["GSEA_r"]]
#GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
#type <- snakemake@wildcards[["msign"]]
#gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
#gene_res_df <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.gseain_2.tsv", 
#                          quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici



em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)

#save.image("gsea_results_c6.R")

#GSEA_r <- write.table(em@result, quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(em@result, file = results, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

#ridgeplot(em, showCategory = 20)
#ggsave(GSEA_ridgeplot, width = 300, height = 107, useDingbats=FALSE, units = "mm")

save.image(imagine)
