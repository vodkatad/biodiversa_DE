## metadata 

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
meta <- snakemake@output[["meta"]]

# meda <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMX_BASALE", type))
meda_f$model <- substr(meda_f$sample_id_R, 1,7)
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))

##removed CRC1502 11/07 post riunione chemio

magnifici_30_resistant <- c("CRC0029", "CRC1078", "CRC0382", "CRC0031", "CRC0772", "CRC0077", "CRC0534",
                            "CRC0151", "CRC0479", "CRC0137", "CRC0204", "CRC0370", "CRC0568", "CRC0019")
magnifici_30_sens <- c("CRC0322", "CRC0076", "CRC0196", "CRC0542", "CRC0330", "CRC0096", "CRC0297", "CRC0729","CRC0161",
                       "CRC0152","CRC0059","CRC0121","CRC0069","CRC0743","CRC0197")

meda_f_res <- meda_f[meda_f$model %in% magnifici_30_resistant,]
meda_f_res$type <- NULL
meda_f_res$type <- "resistant"

meda_f_sens <- meda_f[meda_f$model %in% magnifici_30_sens,]
meda_f_sens$type <- NULL
meda_f_sens$type <- "sensitive"

meda_fin <- rbind(meda_f_res, meda_f_sens)
col_order <- c("model", "batch", "type", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")


write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)
#mmmm <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/metadata_biobanca.tsv"
#mmmmm <- read.table(mmmm, quote = "", sep = "\t", header = TRUE)
