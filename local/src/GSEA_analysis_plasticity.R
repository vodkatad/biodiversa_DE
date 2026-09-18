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

gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)

if (type == "plasticity") {
  sign_f <- '/mnt/trcanmed/snaketree/stash/collection_signature_def.tsv'
  sign <- read.table(sign_f, sep = "\t", header = TRUE)
  sign <- sign[sign$Kept, ]
  m_t2g <- sign[, c('Annotation', 'Gene')]
  colnames(m_t2g) <- c('gs_name', 'gene_symbol')
} else {
  m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>%
    dplyr::select(gs_name, gene_symbol)
}

em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)

write.table(em@result, file = GSEA_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

#ridgeplot(em, showCategory = 20)
res <- em@result
res$sign <- factor(ifelse(res$p.adjust < 0.05, "Yes", "No"), levels = c("Yes", "No"))

showCategory <- 20
res <- res[order(res$p.adjust, decreasing = FALSE), ]
res <- head(res, showCategory)

gs2id <- lapply(res$core_enrichment, function(x) strsplit(x, "/")[[1]])
names(gs2id) <- res$Description

geneList <- em@geneList

gs2val.df <- do.call(rbind, lapply(names(gs2id), function(nm) {
  data.frame(category = nm, value = geneList[gs2id[[nm]]])
}))

gs2val.df$category <- factor(gs2val.df$category, levels = rev(res$Description))

gs2val.df <- merge(gs2val.df, res[, c("Description", "sign")],
                   by.x = "category", by.y = "Description")

label_func <- function(x) stringr::str_wrap(x, width = 30)

ggplot(gs2val.df, aes(x = value, y = category, fill = sign)) +
  ggridges::geom_density_ridges() +
  scale_fill_manual(
    name = "p.adj < 0.05",
    values = c("Yes" = "darkgreen", "No" = "grey70")
  ) +
  scale_y_discrete(labels = label_func) +
  xlab(NULL) + ylab(NULL) +
  DOSE::theme_dose()

ggsave(GSEA_ridgeplot, width = 300, height = 107, useDingbats = FALSE, units = "mm")

save.image('GSEA.Rdata')