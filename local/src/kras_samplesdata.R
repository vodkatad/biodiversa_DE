## metadata 

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
casi <- snakemake@input[["txt"]]
meta <- snakemake@output[["meta"]]

#casi <- "/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/rasdep_20220207.txt"
casi_f <- read.table(casi, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi_f <- casi_f %>% mutate(quartile = ntile(rank, 4))
casi_f <- casi_f %>% filter(quartile == 1 | quartile == 4)
casi_4 <- casi_f$smodel[casi_f$quartile == 4]
casi_1 <- casi_f$smodel[casi_f$quartile == 1]

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

meda_f_4 <- meda_f[meda_f$model %in% casi_4,]
meda_f_4$type <- NULL
meda_f_4$type <- "4Q"

meda_f_1<- meda_f[meda_f$model %in% casi_1,]
meda_f_1$type <- NULL
meda_f_1$type <- "1Q"

meda_fin <- rbind(meda_f_4, meda_f_1)
col_order <- c("model", "batch", "type", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
names(meda_fin)[names(meda_fin) == 'model'] <- 'sample'

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)
