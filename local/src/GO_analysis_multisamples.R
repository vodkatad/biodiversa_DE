library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)


gene_res_f <- snakemake@input[["gene_list"]]     
gene_univ_f <- snakemake@input[["gene_univ"]]    
GO_r <- snakemake@output[["GO_r"]]               
out_dir <- snakemake@output[["out_dir"]]         
direction <- snakemake@wildcards[["direction"]]  

dirs_to_create <- unique(c(dirname(GO_r), out_dir))
for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = FALSE)
gene_univ_df <- read.table(gene_univ_f, quote = "", sep = "\t", header = FALSE)

geneList <- as.character(gene_res_df$V1)
geneUni <- as.character(gene_univ_df$V1)


egocc <- enrichGO(gene = geneList,
                  universe = geneUni,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = FALSE)

egomf <- enrichGO(gene = geneList,
                  universe = geneUni,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = FALSE)

egobp <- enrichGO(gene = geneList,
                  universe = geneUni,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = FALSE)


ggsave(file.path(out_dir, paste0("CC_", direction, ".pdf")),
       plot = barplot(egocc, showCategory = 20),
       width = 300, height = 107, units = "mm", useDingbats = FALSE)

ggsave(file.path(out_dir, paste0("MF_", direction, ".pdf")),
       plot = barplot(egomf, showCategory = 20),
       width = 300, height = 107, units = "mm", useDingbats = FALSE)

ggsave(file.path(out_dir, paste0("BP_", direction, ".pdf")),
       plot = barplot(egobp, showCategory = 20),
       width = 300, height = 107, units = "mm", useDingbats = FALSE)


egocc@result$ontology <- "CC"
egomf@result$ontology <- "MF"
egobp@result$ontology <- "BP"

egoall_df <- rbind(egocc@result, egomf@result, egobp@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method = "BH")

write.table(egoall_df,
            file = GO_r,
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

save.image(file.path(out_dir, paste0("GO_", direction, ".Rdata")))