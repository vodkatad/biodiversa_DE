#print(snakemake@params[["threads"]]) # rid to name
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(showtext)
library(ggrastr)
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- (size) * death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize),#, family='Arial'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size=0.508/0.564), # origin of this ratio is honestly not known, empirical
  axis.ticks = element_line(color = "black", size=0.508/0.564),
  axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),
  panel.background = element_blank()
)
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

#geni <- c("LCN2","DMBT1","ITLN1","NOS2","SERPINE1","VIM","COL6A1","COL7A1","LAMA5", "STAG3", "SPINK4", "FCGBP", "AL354836.1")
geni <- c("VIM", "COL7A1", "COL6A1", "SERPINE1", "LAMA5", "LCN2", "DMBT1", "ITLN1", "NOS2", "CAPN9", "SPINK4", "FCGBP", "MFAP4", "VGLL3", "SPON2", "OLFML2A", "EMP3", "STAG3")
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
    rasterise(geom_point(aes(col = sign),size=1), dpi=300) + theme_bw() +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE) + # red orange green black -> orange blue green gray 
    ggtitle(title)+unmute_theme
  
  # nsign <- nrow(resnona[resnona$sig=="both",])
  # if (nsign > 20) {
  #   p + geom_text_repel(data=resnona[1:10,], aes(label=rownames(resnona)[1:10]))
  # } else {
  #   p + geom_text_repel(data=resnona[resnona$sig=="both",], aes(label=rownames(resnona[resnona$sig=="both",])))
  # }
  resnona$scelti <- ifelse(rownames(resnona)%in%geni, "YES", "NO")
  resnona$scelti <- factor(resnona$scelti, levels=c("YES", "NO"))
  
  selected_points <- resnona[resnona$scelti == "YES",]

# Add arrows from points to labels
  p + geom_segment(data = selected_points, aes(x = log2FoldChange, y = -log10(padj), 
                   xend = log2FoldChange, yend = -log10(padj) + 1),  # Adjust label position
               color = "black", linetype = "dashed") +  # Dashed lines
  geom_text_repel(data = selected_points, aes(label = rownames(selected_points)), nudge_y = 1.5)
  #p + geom_text_repel(data=resnona[resnona$scelti=="YES",], aes(label=rownames(resnona[resnona$scelti=="YES",])))
  ggsave(outfile, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")
}


p <- plot_volcano(resnona_df, alpha, lfc, volcano, title)


write.table(resnona_df, file=tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

save.image(paste0(tsv, "_DESeq.Rdata"))
