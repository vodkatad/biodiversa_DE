library(tidyverse)


metadata_o_f <- snakemake@input[["metadata"]] 
samples <- snakemake@input[["original_w3"]]
meta <- snakemake@output[["sample"]]

#samples <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
s <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
s$batch <- NULL
s <- s[!duplicated(s$sample),]

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMO_BASALE", type))
meda_f$sample <- substr(meda_f$sample_id_R, 1,7)
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_f$type <- NULL

merged <- merge(s, meda_f, by = "sample")
res <- merged
res <- res[!duplicated(res$sample_id_R),]
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)