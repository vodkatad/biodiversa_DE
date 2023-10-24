library(tidyverse)

deg_f <- snakemake@input[["data"]]
geni_f <- snakemake@input[["genes_data"]]
results <- snakemake@output[["tsv"]]
results_filter <- snakemake@output[["tsv_filter"]]

#geni <- read.table("/scratch/trcanmed/DE_RNASeq/local/share/data/HRgenes_all.txt",quote = "", sep = "\t", header = F, stringsAsFactors = F)

geni <- read.table(geni_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#deg <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv",quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg <- read.table(deg_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#save.image("muto.Rdata")

geni$geni <- paste0("H_", geni$V3)
geni <- geni$geni

deg$genes <- rownames(deg)

deg <- deg %>% filter(genes %in% geni)

deg$correction_padj <- p.adjust(deg$pvalue, method = 'fdr' )
deg$genes <- NULL

genes <- rownames(deg)
rownames(deg) <- NULL
deg <- cbind(genes,deg)
deg <- deg %>% mutate(genes = gsub("H_", "", genes))

write.table(deg, file=results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

deg2 <- deg %>% filter(pvalue < 0.05)
deg2 <- deg2[,c(1,3,6)]

write.table(deg2, file=results_filter, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
#v <- gsub("H_", "", rownames(deg))
#write.table(v, "/scratch/trcanmed/DE_RNASeq/local/share/data/livio2", quote = FALSE, sep = "\t",
#            row.names = FALSE, col.names = FALSE)