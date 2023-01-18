### samples data for DEG tcf7l2
#4 tcl knock out resistant
#1 tcl know out sensibili

library(tidyverse)
library(readxl)

metadata_o_f <- snakemake@input[["metadata"]] 
casi <- snakemake@input[["txt"]]
meta <- snakemake@output[["meta"]]

mut <- read_excel(casi)
names(mut)[names(mut) == "TCF7L2 mut status"] <- "type"
mut$case <- gsub(" ", "", mut$case)

mut <- as.data.frame(mut[,c(1,5)])
for (i in seq(mut$case)) {
  if (mut[i,"type"] == "wt") {
    mut[i, "type"] <- "wt"
  } else {
    mut[i, "type"] <- "mut"
  }
}

names(mut)[names(mut) == "case"] <- "model"

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

meda_fin <- merge(mut, meda_f, by = "model")

col_order <- c("model", "batch", "type.x", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
names(meda_fin)[names(meda_fin) == 'model'] <- 'sample'
names(meda_fin)[names(meda_fin) == "type.x"] <- "type"

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)