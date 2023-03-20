## samples datas for LMx vs LMO per interferoni marika

library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
d <- snakemake@input[["txt"]]
meta <- snakemake@output[["meta"]]

#d <- "/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/samples_data"
d <- read.table(d, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
d <- d %>% filter(type == "Dep_1Q")
sample <- d$sample

#meta_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda <- meda %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda$ShortID <- substr(meda$sample_id_R, 1, 7)
meda_x <- filter(meda, grepl("LMX_BASALE", type))
meda_o <- filter(meda, grepl("LMO_BASALE", type))

meda_x <- meda_x %>% filter(ShortID %in% sample)
meda_o <- meda_o %>% filter(ShortID %in% sample)

meda_x <- meda_x[,c(1, 6, 7, 10)]
meda_o <- meda_o[,c(1, 6, 7, 10)]

meda_x$type <- "LMX"
meda_o$type <- "LMO"

result <- as.data.frame(rbind(meda_o, meda_x))
rownames(result) <- result$sample_id_R
result$sample_id_R <- NULL

write.table(result, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)