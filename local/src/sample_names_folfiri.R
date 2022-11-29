## deg livio per treated models

library(tidyverse)


metadata_o_f <- snakemake@input[["metadata"]] 
casi_f <- snakemake@input[["folfiri"]]
meta <- snakemake@output[["sample"]]


#casi_f <- "/scratch/trcanmed/pdxopedia/local/share/data/treats/september2022/folfiri-cetuxi_3w_treated_models.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi <- casi[order(casi$folfiri.vol_3w, decreasing = TRUE),]
casi$ntile <- ntile(casi$folfiri.vol_3w, 4)
casi$ntile3 <- ntile(casi$folfiri.vol_3w, 3)

ggplot(casi, aes(x = reorder(model, -folfiri.vol_3w), y = folfiri.vol_3w, fill = as.factor(ntile))) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, size = 5))

ggplot(casi, aes(x = reorder(model, -folfiri.vol_3w), y = folfiri.vol_3w, fill = ntile3)) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, size = 5))

quarti <- quantile(casi$folfiri.vol_3w, c(1/3, 2/3))

casi$quartile <- NA

casi$quartile <- ifelse(casi$folfiri.vol_3w < quarti[1], 1, ifelse(casi$folfiri.vol_3w > quarti[2], 3, 2))                

casi <- casi %>% filter(quartile == 1 | quartile == 3)

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
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
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))

merged <- merge(casi, meda_f, by = "model")
res <- as.data.frame(merged[, c(1, 10, 11, 12)])
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
names(res)[names(res) == "model"] <- "sample"

for (i in seq(rownames(res))){
  if (res[i, "quartile"] == 1) {
    res[i, "quartile"] <- "PR"
  } else {
    res[i, "quartile"] <- "PD"
  }
}

names(res)[names(res) == "quartile"] <- "type"

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)