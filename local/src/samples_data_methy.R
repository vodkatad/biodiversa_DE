library(tidyverse)

meta_f <- snakemake@input[["metadata"]]
methy <- snakemake@input[["metilation"]]
results <- snakemake@output[["sample"]]

#meta_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meta_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda <- meda %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda$model <- substr(meda$sample_id_R, 1, 7)
meda <- filter(meda, grepl("LMX_BASALE", type))

#methy <- "/scratch/trcanmed/pdx_methylation/local/share/data/annotations/All_samples_info_final_060422_wclusters-cn_CRISvsd.tsv"
m <- read.table(methy, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
m$class <- substr(m$type, 8, 10) 
m <- m %>% filter(class == "LMX")

m2 <- as.data.frame(cbind(m$model, m$cluster))
colnames(m2) <- c("model", "cluster")

merged <- merge(meda, m2, by = "model")

res <- as.data.frame(cbind(merged$model, merged$sample_id_R, merged$batch, merged$cluster))
colnames(res) <- c("model", "sample_id_R", "batch", "cluster_n")
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL

res$cluster <- NA
for (i in seq(length(res$model))){
  if (res[i, "cluster_n"] == 1) {
    res[i, "cluster"] <- "clu_1"
  } else if (res[i, "cluster_n"] == 2) {
    res[i, "cluster"] <- "clu_2"
  } else if (res[i, "cluster_n"] == 3) {
    res[i, "cluster"] <- "clu_3"
  } else if (res[i, "cluster_n"] == 4) {
    res[i, "cluster"] <- "clu_4"
  } else {
    res[i, "cluster"] <- "clu_5"
  }
}

res$cluster_n <- NULL

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

