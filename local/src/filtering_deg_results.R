## heatmap cluster methy

library(tidyverse)
library(pheatmap)

res <- snakemake@input[["deg"]]
vsd <- snakemake@input[["expr"]]
data <- snakemake@input[["sample"]]
table1 <- snakemake@output[["sign"]]
table2 <- snakemake@output[["top"]]
heat1 <- snakemake@output[["pbig"]]
heat2 <- snakemake@output[["psmall"]]

#res <- "/scratch/trcanmed/DE_RNASeq/dataset/methy_clusters/cluster_methy_cutoff0.05.deseq2.tsv"
res <- read.table(res, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

fil <- res
fil <- fil %>% filter(padj < 0.001)
fil$genes <- rownames(fil)

#vsd <- "/scratch/trcanmed/DE_RNASeq/dataset/methy_clusters/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$genes <- rownames(vsd)

merged <- merge(vsd, fil, by = "genes")
rownames(merged) <- merged$genes
merged$genes <- NULL
merged$log2FoldChange <- NULL
merged$baseMean <- NULL
merged$lfcSE <- NULL
merged$padj <- NULL
merged$pvalue <- NULL
merged$stat <- NULL
tvsd <- as.data.frame(t(merged))
#tvsd$model <- substr(rownames(tvsd), 1, 7)
tvsd$genealogy <- rownames(tvsd)

#data <- "/scratch/trcanmed/DE_RNASeq/dataset/methy_clusters/samples_datas"
data <- read.table(data, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
data$batch <- NULL
data$genealogy <- rownames(data)

merged2 <- merge(tvsd, data, by = "genealogy")
rownames(merged2) <- merged2$genealogy
merged2$genealogy <- NULL
merged2$model <- NULL

save.image("methy.Rdata")
### fare la media per gene dividendo per cluster

cluster_pigri <- function(clu, df) {
  df <- df %>% filter(cluster == clu)
  df$cluster <- NULL
  df <- as.data.frame(t(df), stringsAsFactors = FALSE)
  name <- paste0(clu, "_means")
  df[, name] <- rowMeans(df)
  df$genes <- rownames(df)
  df <- select(df, name, genes)
  return(df)
}
#df <- merged2
#clu <- "clu_1"
clu1 <- cluster_pigri("clu_1", merged2)
clu2 <- cluster_pigri("clu_2", merged2)
clu3 <- cluster_pigri("clu_3", merged2)
clu4 <- cluster_pigri("clu_4", merged2)
clu5 <- cluster_pigri("clu_5", merged2)

clu12 <- merge(clu1, clu2, by = "genes")
clu123 <- merge(clu12, clu3, by = "genes")
clu1234 <- merge(clu123, clu4, by = "genes")
clu_fin <- merge(clu1234, clu5, by = "genes")
rownames(clu_fin) <- clu_fin$genes
clu_fin$genes <- NULL

## mettere qui una write.table
write.table(clu_fin, file = table1, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

png(heat1)
pheatmap(clu_fin, cluster_rows = FALSE, cluster_cols = FALSE, labels_row = "cases")
dev.off()
#riga per riga calcolare la ds/media i valori più altri hanno più facilità a variare prendo poi i top 100 ordinati per questa media

clusters <- clu_fin
clusters$means <- rowMeans(clusters)
clusters$sd <- apply(clusters, 1, sd)   
clusters$formula <- (clusters$sd/clusters$means)
clusters <- clusters[order(-clusters$formula),]
res <- head(clusters, 100)
results <- res[1:5]

write.table(results, file = table2, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

png(heat2)
pheatmap(results, cluster_rows = FALSE, cluster_cols = FALSE, labels_row = "cases")
dev.off()