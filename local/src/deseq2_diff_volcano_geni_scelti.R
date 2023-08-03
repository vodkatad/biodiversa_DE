#print(snakemake@params[["threads"]]) # rid to name
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)

register(MulticoreParam(as.numeric(snakemake@params[["threads"]])))
threads <- as.numeric(snakemake@params[["threads"]])
parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    register(MulticoreParam(threads))
    parallel <- TRUE
}
alpha <- as.numeric(snakemake@params[["alpha"]])
lfc <- as.numeric(snakemake@params[["lfc"]]) # used only for volcano plots, the tsv printed lists all non NA results!

#print(snakemake@input[[1]]) # .RData
#print(snakemake@params[["class"]]) # which columns need to be compared
#print(snakemake@params[["nom"]]) # this vs
#print(snakemake@params[["den"]]) # this one
con <- c(snakemake@params[["factor"]], snakemake@params[["nom"]], snakemake@params[["den"]])
volcano <- snakemake@output[["volcano"]]
tsv <- snakemake@output[["tsv"]]
#load overwrites our snakemake object thus we need to put aside our parameters before.

save.image(paste0(tsv, "_DESeq.Rdata"))

load(snakemake@input[[1]])

res <- results(dds, alpha=alpha, contrast=con, parallel=parallel)

resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])
title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])

geni <- c("LCN2","DMBT1","ITLN1","NOS2","SERPINE1","VIM","COL6A1","COL7A1","LAMA5")
genes_or <- rownames(resnona_df)
rownames(resnona_df) <- NULL
resnona_df <- cbind(genes_or,resnona_df)
resnona_df <- resnona_df %>% mutate(genes_or = gsub("H_", "", genes_or))
rownames(resnona_df) <- resnona_df$genes_or
resnona_df$genes_or <- NULL

#prova <- resnona_df
plot_volcano <- function(resnona, alpha, lfc, outfile, title) {
  #title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
  #resnona <- res[!is.na(res$padj),]
  resnona$sign <- ifelse(abs(resnona$log2FoldChange) > lfc & resnona$padj < alpha, "both", ifelse(abs(resnona$log2FoldChange) > lfc, "LFC",
                                                                                                  ifelse(resnona$padj < alpha, "padj", "NS")))
  resnona$sign <- factor(resnona$sig, levels=c("LFC", "padj", "both", "NS"))
  resnona[resnona$padj ==0,"padj"] <- .Machine$double.xmin
  p <- ggplot(resnona, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(col = sign),size=0.5) + theme_bw() +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE) + # red orange green black -> orange blue green gray 
    ggtitle(title)
  
  # nsign <- nrow(resnona[resnona$sig=="both",])
  # if (nsign > 20) {
  #   p + geom_text_repel(data=resnona[1:10,], aes(label=rownames(resnona)[1:10]))
  # } else {
  #   p + geom_text_repel(data=resnona[resnona$sig=="both",], aes(label=rownames(resnona[resnona$sig=="both",])))
  # }
  resnona$scelti <- ifelse(rownames(resnona)%in%geni, "YES", "NO")
  resnona$scelti <- factor(resnona$scelti, levels=c("YES", "NO"))
  p + geom_text_repel(data=resnona[resnona$scelti=="YES",], aes(label=rownames(resnona[resnona$scelti=="YES",])))
  ggsave(outfile)
}


plot_volcano(resnona_df, alpha, lfc, volcano, title)


write.table(resnona_df, file=tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

save.image(paste0(tsv, "_DESeq.Rdata"))
