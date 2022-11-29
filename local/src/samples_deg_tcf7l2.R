### samples data for DEG tcf7l2
#4 tcl knock out resistant
#1 tcl know out sensibili

library(tidyverse)
library(readxl)

metadata_o_f <- snakemake@input[["metadata"]] 
casi <- snakemake@input[["txt"]]
meta <- snakemake@output[["meta"]]

mut <- read_excel(casi)
names(mut)[names(mut) == "Co-co score N2 sgRNA"] <- "co_comp_N2"
mut$case <- gsub(" ", "", mut$case)

mut <- as.data.frame(mut[,c(1,3)])
mut <- mut[order(mut[,2], decreasing = TRUE),]
quarti <- quantile(mut$co_comp_N2, c(1/4, 3/4))
mut$quartile <- NA
rownames(mut) <- mut$case
mut$quartile <- ifelse(mut$co_comp_N2 < quarti[1], 1, ifelse(mut$co_comp_N2 > quarti[2], 4, 23))

#write.xlsx(mut, file = "tcf7l2_order_quarti.xlsx")

mut <- mut %>% filter(quartile == 1| quartile == 4)
casi_res4 <- mut$case[mut$quartile == 4]
casi_sen1 <- mut$case[mut$quartile == 1]

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

meda_f_4 <- meda_f[meda_f$model %in% casi_res4,]
meda_f_4$type <- NULL
meda_f_4$type <- "Res_4Q"

meda_f_1<- meda_f[meda_f$model %in% casi_sen1,]
meda_f_1$type <- NULL
meda_f_1$type <- "Sen_1Q"

meda_fin <- rbind(meda_f_4, meda_f_1)
col_order <- c("model", "batch", "type", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
names(meda_fin)[names(meda_fin) == 'model'] <- 'sample'

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)