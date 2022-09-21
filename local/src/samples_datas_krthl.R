### samples data deg krt high e low

library(tidyverse)

krt_f <- snakemake@input[["class_krt"]]
meda <- snakemake@input[["metadata"]]
results <- snakemake@output[["sample"]]

#krt_f <- "/home/mferri/Pd_high_low_krt.tsv"
krt <- read.table(krt_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
krt$sample_id_R <- rownames(krt)

#meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_f <- meda_f %>% filter(sample_id_R %in% krt$sample_id_R)

merged <- merge(krt, meda_f, by = "sample_id_R")
res <- as.data.frame(cbind(merged$sample_id_R, merged$batch, merged$KRT_classification))
colnames(res) <- c("genealogy", "batch", "class")
res$model <- substr(res$genealogy, 0, 7)
rownames(res) <- res$genealogy
res$genealogy <- NULL
res <- res[,c(3,1,2)]
res$model <- as.character(res$model)
res$batch <- as.factor(res$batch)
res$class <- as.character(res$class)

for (i in rownames(res)) {
  if (res[i, "class"] == "KRT_high") {
    res[i, "class"] <- "H"
  } else if (res[i, "class"] == "KRT_low") {
    res[i, "class"] <- "L"
  } else {
    res[i, "class"] <- "I"
  }
}

write.table(res, results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)