## heatmap fpkm log+1 

library(tidyverse)
library(pheatmap)

fpkm_f <- snakemake@input[["data"]]
samples_f <- snakemake@input[["samples_data"]]
results <- snakemake@output[["heatmap"]]

#fpkm_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/fpkm.tsv.gz"
fpkm <- read.table(fpkm_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(fpkm)
rownames(fpkm) <- NULL
fpkm <- cbind(genes,fpkm)
fpkm <- fpkm %>% mutate(genes = gsub("H_", "", genes))

geni_livio <- c("FBXO18", "SUMO1", "UBE2I", "BRCA1", "RRP1",
                "FBXO5", "RING1", "RFWD3",
                "UCHL3", "PARPBP", "BLM")

fpkm_geni <- fpkm %>% filter(genes %in% geni_livio)                
rownames(fpkm_geni) <- fpkm_geni$genes
fpkm_geni$genes <- NULL

fpkm_log <- log(fpkm_geni+1) 
fpkm_log <- as.data.frame(t(fpkm_log))

#samples_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
samples <- read.table(samples_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- samples[order(samples$type),]
samples$sample <- NULL
samples$batch <- NULL

samples$genealogy <- rownames(samples)
fpkm_log$genealogy <- rownames(fpkm_log)
fpkm_log <- merge(fpkm_log, samples, by="genealogy")
fpkm_log <- fpkm_log[order(fpkm_log$type),]
rownames(fpkm_log) <- fpkm_log$genealogy
fpkm_log$genealogy <- NULL
fpkm_log$type <- NULL
samples$genealogy <- NULL
fpkm_log <- fpkm_log[, geni_livio]

an_col <- as.data.frame(matrix(ncol = 1, nrow = 11))
rownames(an_col) <- colnames(fpkm_log)
an_col$V1 <- c(rep("Ubiquitin ligase",8), "Deubiquitinase", rep("Elicase",2))
names(an_col)[names(an_col)=="V1"] <- "Enzyme"

### scale fa per ogni colonna elemento - media/ sd
p <- pheatmap(fpkm_log, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = samples, annotation_col = an_col, scale = "column")
ggsave(p, filename = results)
