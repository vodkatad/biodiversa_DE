## connector Aa vs all 2nd round post revision 

library(tidyverse)

connector_f <- snakemake@input[["connector_class"]]
meta_f <- snakemake@input[["metadata"]]
uno <- snakemake@wildcards[['first']]
due <- snakemake@wildcards[['second']]
results <- snakemake@output[["sample"]]

#connector_f <- "/scratch/trcanmed/connector/local/share/data/Models.RDs"
connector <- readRDS(connector_f)
#connector$group <- substr(connector$col, 1,1)
connector <- connector %>% filter(Cluster %in% c(uno, due))
names(connector)[names(connector) == "ShortID"] <- "model"

#meta_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meta_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda <- meda %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda$model <- substr(meda$sample_id_R, 1, 7)
meda <- filter(meda, grepl("LMX_BASALE", type))

merged <- merge(connector, meda, by = "model")
merged <- merged[, c("sample_id_R", "model", "Cluster", "batch")]
rownames(merged) <- merged$sample_id_R
merged$sample_id_R <- NULL

write.table(merged, results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)