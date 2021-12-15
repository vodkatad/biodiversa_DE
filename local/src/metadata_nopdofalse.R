## metadata 

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]]
pdo_buoni_or <- snakemake@input[["pdo"]]
meta <- snakemake@output[["meta"]]

#meda <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#pdo_buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
print(pdo_buoni_or)
pdo_buoni <- read.table(pdo_buoni_or, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo_buoni_f <- filter(pdo_buoni, buoni == FALSE)
colnames(pdo_buoni_f) <- c("sample", "buoni")

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
meda_xh <- meda_f
meda_xh <- filter(meda_xh, type == "LMX_BASALE" | type == "LMH")
meda_xh <- meda_xh %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_xh <- meda_xh %>% remove_rownames %>% column_to_rownames(var="sample_id_R")

## rimozione 3 pdo malefici
meda_o <- meda_f
meda_o <- filter(meda_o, type == "LMO_BASALE")
meda_o <- meda_o %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_o <- meda_o %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
mergedo_good <- merge(meda_o, pdo_buoni_f, by= "sample")
removal_sample <- mergedo_good$sample

meda_o <- meda_o[!meda_o$sample %in% removal_sample,]

## rimozione CRC0177LMO0D..
meda_o <- tibble::rownames_to_column(meda_o, "genealogy")
crc177 <- c("CRC0177LMO0D04021002R01000")
meda_o <- meda_o[!meda_o$genealogy %in% crc177,]
meda_o <- meda_o %>% remove_rownames %>% column_to_rownames(var="genealogy")

meda_fin <- rbind(meda_xh, meda_o)


write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)

