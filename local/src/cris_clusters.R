### DEG with cris classes 

library(tidyverse)

meta_f <- snakemake@input[["metadata"]]
cris_f <- snakemake@input[["cris"]]
results <- snakemake@output[["sample"]]


#meta_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meta_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda <- meda %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda$model <- substr(meda$sample_id_R, 1, 7)
meda <- filter(meda, grepl("LMX_BASALE", type))

#cris_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/model_cris-right.tsv"
cris <- read.table(cris_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cris$CRIS_PDO <- NULL
cris <- cris %>% filter(!CRIS_PDX == "NC")
cris <- cris %>% mutate(CRIS_PDX = gsub("-", "_", CRIS_PDX))

merged <- merge(meda, cris, by = "model")

res <- as.data.frame(cbind(merged$model, merged$sample_id_R, merged$batch, merged$CRIS_PDX))
colnames(res) <- c("model", "genealogy", "batch", "CRIS_PDX")
rownames(res) <- res$genealogy
res$genealogy <- NULL

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)