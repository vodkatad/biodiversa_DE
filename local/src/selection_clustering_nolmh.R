library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]]
pdo_buoni_or <- snakemake@input[["pdo"]]
meta <- snakemake@output[["meta"]]

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#pdo_buoni_or <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo_buoni <- read.table(pdo_buoni_or, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo_buoni_f <- filter(pdo_buoni, buoni == FALSE)
#colnames(pdo_buoni_f) <- c("sample", "buoni")
pdo_buoni_t <- filter(pdo_buoni, buoni == TRUE)
colnames(pdo_buoni_t) <- c ("sample", "buoni")

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

## tengo solo i PDO true così sono certa siano i PDO validati
## (eliminando solo i falsi tengo più PDO di quelli effettivamente validati)
meda_o <- meda_f
meda_o <- filter(meda_o, type == "LMO_BASALE")
meda_o <- meda_o %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_o <- meda_o %>% remove_rownames %>% column_to_rownames(var="sample_id_R")

mergedo_good <- merge(meda_o, pdo_buoni_t, by= "sample")
keep_sample_o <- mergedo_good$sample

meda_o <- meda_o[meda_o$sample %in% keep_sample_o,]

## rimozione CRC0177LMO0D..
meda_o <- tibble::rownames_to_column(meda_o, "genealogy")
crc177 <- c("CRC0177LMO0D04021002R01000")
meda_o <- meda_o[!meda_o$genealogy %in% crc177,]
meda_o <- meda_o %>% remove_rownames %>% column_to_rownames(var="genealogy")

## vedere quanti LMX elimino vedendo la corrispondenza con gli organoidi buoni
meda_x <- meda_f
meda_x <- filter(meda_x, type == "LMX_BASALE")
meda_x <- meda_x %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_x <- meda_x %>% remove_rownames %>% column_to_rownames(var="sample_id_R")

mergedx_good <- merge(meda_x, pdo_buoni_t, by= "sample")
keep_sample_x <- mergedx_good$sample

meda_x <- meda_x[meda_x$sample %in% keep_sample_x,]

o <- meda_o$sample
x <- meda_x$sample

difference <- setdiff(o,x)

meda_fin <- rbind(meda_o, meda_x)
meda_fin <- meda_fin[!meda_fin$sample %in% difference,]

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)