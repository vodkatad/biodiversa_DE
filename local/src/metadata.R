## metadata 

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
meta <- snakemake@output[["meta"]]

#meda <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL


meda_f$sample <- substr(meda_f$sample_id_R, 1, 7)
col_order <- c("sample", "batch", "type", "sample_id_R")
meda_f <- meda_f[, col_order]
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))
meda_f <- filter(meda_f, type == "LMO_BASALE" | type == "LMX_BASALE" | type == "LMH")
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_f <- meda_f %>% remove_rownames %>% column_to_rownames(var="sample_id_R")

write.table(meda_f, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)

