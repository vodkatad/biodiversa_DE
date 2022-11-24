## metadata 

library(tidyverse)
library(readxl)

metadata_o_f <- snakemake@input[["metadata"]] 
casi <- snakemake@input[["txt"]]
meta <- snakemake@output[["meta"]]

#casi <- "/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/rasdep_20220207.txt"
#casi_f <- read.table(casi, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi_f <- as.data.frame(read_excel(casi))
casi_f$...3 <- NULL
casi_f$...4 <- NULL
casi_f$...5 <- NULL
casi_f$...6 <- NULL
casi_f$...7 <- NULL
quarti <- quantile(casi_f$RASDEPscore, c(1/4, 3/4))
casi_f$quartile <- NA
rownames(casi_f) <- casi_f$Case
#casi_f$RASDEPscore <- as.integer(casi_f$RASDEPscore)

casi_f$quartile <- ifelse(casi_f$RASDEPscore < quarti[1], 1, ifelse(casi_f$RASDEPscore > quarti[2], 4, 23))

# for (i in rownames(casi_f)) {
#   if(casi_f[i, "RASDEPscore"] == 1.065944 | casi_f[i, "RASDEPscore"] > 0.4704496){
#     casi_f[i, "quartile"] <- 4
#   } else if (casi_f[i, "RASDEPscore"] < 0.4704496 & casi_f[i, "RASDEPscore"] > 0.2570949){
#     casi_f[i, "quartile"] <- 3
#   } else if (casi_f[i, "RASDEPscore"] < 0.2570949 & casi_f[i, "RASDEPscore"] > 0.08588025){
#     casi_f[i, "quartile"] <- 2
#   } else {
#     casi_f[i, "quartile"] <- 1
#   }
# }

#casi_f <- casi_f %>% mutate(quartile = ntile(RASDEPscore, 4))
casi_f <- casi_f %>% filter(quartile == 1| quartile == 4)
casi_ind4 <- casi_f$Case[casi_f$quartile == 4]
casi_dip1 <- casi_f$Case[casi_f$quartile == 1]

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

meda_f_4 <- meda_f[meda_f$model %in% casi_ind4,]
meda_f_4$type <- NULL
meda_f_4$type <- "Ind_4Q"

meda_f_1<- meda_f[meda_f$model %in% casi_dip1,]
meda_f_1$type <- NULL
meda_f_1$type <- "Dep_1Q"

meda_fin <- rbind(meda_f_4, meda_f_1)
col_order <- c("model", "batch", "type", "sample_id_R")
meda_fin <- meda_fin[, col_order]
meda_fin <- meda_fin %>% remove_rownames %>% column_to_rownames(var="sample_id_R")
names(meda_fin)[names(meda_fin) == 'model'] <- 'sample'

write.table(meda_fin, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)
