## metadata 

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
meta <- snakemake@output[["meta"]]

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMO_BASALE", type))
meda_f$model <- substr(meda_f$sample_id_R, 1,7)
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))

palbo_res <- c("CRC1257", "CRC1588", "CRC1589")
palbo_non_res <- c("CRC0078", "CRC0322", "CRC0534")

meda_f_res <- meda_f[meda_f$model %in% palbo_res,]
meda_f_res$type <- NULL
meda_f_res$type <- "responder"

meda_f_non_res <- meda_f[meda_f$model %in% palbo_non_res,]
meda_f_non_res$type <- NULL
meda_f_non_res$type <- "non_responder"

meda_fin <- rbind(meda_f_res, meda_f_non_res)
col_order <- c("model", "batch", "type", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
names(meda_fin)[names(meda_fin) == 'model'] <- 'sample'

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)
