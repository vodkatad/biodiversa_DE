library(tidyverse)

recidive <- snakemake@input[["recidive_f"]]
metadata_o_f <- snakemake@input[["metadata"]] 
meta <- snakemake@output[["sample"]]

#recidive <- "/home/mferri/recidive_classificazione.tsv"
recidive <- read.table(recidive, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
recidive <- recidive[,c("codice_candiolo", "recidiva_class")]
recidive <- recidive %>% filter(!recidiva_class == "never")
rownames(recidive) <- recidive$codice_candiolo
for (i in rownames(recidive)) {
  if (recidive[i, "recidiva_class"] == "<=6") {
    recidive[i, "type"] <- "precoce"
  } else {
    recidive[i, "type"] <- "noprecoce"
  }
}

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL
lmx <- meda_f %>% filter(type %in% c("LMX_BASALE", "LMX_BASALE.1"))
lmx$codice_candiolo <- substr(lmx$sample_id_R, 1, 7)

res <- merge(recidive, lmx, by="codice_candiolo")
res$recidiva_class <- res$type.x
res$type.x <- NULL
res$type.y <- NULL
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)

deg <- "/scratch/trcanmed/DE_RNASeq/dataset/recidive_deg/class_cutoff0.05-precoce.vs.noprecoce.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg <- deg %>% filter(padj < 0.05)
