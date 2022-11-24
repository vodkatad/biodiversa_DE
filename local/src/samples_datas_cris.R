### deg cris-b contro no b

library(tidyverse)

cris_f <- snakemake@input[["cris_class"]]
meda <- snakemake@input[["metadata"]]
results <- snakemake@output[["sample"]]

#cris <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/cris_vsd_ok_prediction_result.tsv"
cris <- read.table(cris_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cases <- cris$sample.names

#meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_f <- meda_f %>% filter(sample_id_R %in% cases)

merged <- merge(cris, meda_f, by.x = "sample.names", by.y = "sample_id_R")
merged <- as.data.frame(cbind(merged$sample.names, merged$batch, merged$predict.label2))
colnames(merged) <- c("genealogy", "batch", "CRIS")
merged$genealogy <- as.character(merged$genealogy)
merged$batch <- as.factor(merged$batch)
merged$CRIS <- as.character(merged$CRIS)
merged$CRIS_def <- ""

for (i in seq(merged$genealogy)) {
  if (merged[i, "CRIS"] == "CRIS-B") {
    merged[i, "CRIS_def"] <- "B"
  } else {
    merged[i, "CRIS_def"] <- "Other"
  }
}

merged$CRIS <- NULL
names(merged)[names(merged) == 'CRIS_def'] <- 'cris'
merged$model <- substr(merged$genealogy, 0 , 7)
rownames(merged) <- merged$genealogy
merged$genealogy <- NULL
merged <- merged[,c(3,1,2)]

write.table(merged, results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)