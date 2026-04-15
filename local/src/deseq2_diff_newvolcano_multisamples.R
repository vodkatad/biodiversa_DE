library(DESeq2)
library(ggplot2)
library(ggrepel)

threads <- as.numeric(snakemake@params[["threads"]])
parallel <- FALSE
if (threads > 1) {
  library("BiocParallel")
  register(MulticoreParam(threads))
  parallel <- TRUE
}

alpha <- as.numeric(snakemake@params[["alpha"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
con <- c(snakemake@params[["factor"]], snakemake@params[["nom"]], snakemake@params[["den"]])
volcano <- snakemake@output[["volcano"]]
tsv <- snakemake@output[["tsv"]]

save.image(paste0(tsv, "_DESeq.Rdata"))

load(snakemake@input[[1]])

res <- results(dds, alpha=alpha, contrast=con, parallel=parallel)

plot_volcano <- function(resnona, alpha, lfc, outfile, title) {
  resnona$sign <- ifelse(resnona$log2FoldChange > lfc & resnona$padj < alpha, "up",
                  ifelse(resnona$log2FoldChange < -lfc & resnona$padj < alpha, "down", "no_diff"))
  
  resnona$signtot <- ifelse(abs(resnona$log2FoldChange) > lfc & resnona$padj < alpha, "both", 
                     ifelse(abs(resnona$log2FoldChange) > lfc, "LFC",
                     ifelse(resnona$padj < alpha, "padj", "NS")))
  
  resnona$padj_capped <- pmax(resnona$padj, 1e-300)
  max_y_val <- max(-log10(resnona$padj_capped), na.rm = TRUE)
  y_limit_upper <- max_y_val * 1.3
  
    p <- ggplot(resnona, aes(log2FoldChange, -log10(padj_capped))) +
        geom_point(aes(col = sign), size=0.7, alpha=0.6) + 
        theme_bw() +
        scale_color_manual(values = c("down" = "blue", "no_diff" = "#999999", "up" = "red"), drop=FALSE) + 
        ggtitle(title) +
        scale_y_continuous(limits = c(0, y_limit_upper)) +
        coord_cartesian(clip = "off") +
        labs(x = "log2FoldChange", y = "-log10(padj)")
  
  nsign <- sum(resnona$signtot == "both")
  repel_params <- list(
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.5,
    force = 2,
    nudge_y = max_y_val * 0.05
  )
  
  label_data <- if (nsign > 20) resnona[1:10,] else resnona[resnona$signtot == "both",]
  p_final <- p + do.call(geom_text_repel, c(list(data=label_data, aes(label=rownames(label_data))), repel_params))
  
  ggsave(outfile, plot = p_final, width = 8, height = 10)
}

resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])
title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
plot_volcano(resnona_df, alpha, lfc, volcano, title)

write.table(resnona_df, file=tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

save.image(paste0(tsv, "_DESeq.Rdata"))