#samples data Ac vs all with recist classification only PD

# <-50% OR
# -%50 35 SD
# +35% PD 

library(tidyverse)

connector_f <- snakemake@input[["connector_class"]]
meta_f <- snakemake@input[["metadata"]]
recist_f <- snakemake@input[["recist_perc_3"]]
uno <- snakemake@wildcards[['first']]
due <- snakemake@wildcards[['second']]
results <- snakemake@output[["sample"]]

#recist_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3.tsv"
recist <- read.table(recist_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
recist$classification <- NA

for (i in seq(length(recist$case))) {
  if (recist[i, "perc"] < -50.0){
    recist[i, "classification"] <- "OR"
  } else if (recist[i, "perc"] > -50.0 & recist[i, "perc"] < 35.0) {
    recist[i, "classification"] <- "SD"
  } else {
    recist[i, "classification"] <- "PD"
  }
}

names(recist)[names(recist) == "case"] <- "ShortID"
recist <- recist %>% filter(classification == "PD")

#connector_f <- "/scratch/trcanmed/connector/local/share/data/tsne_cetuxi_3w.tsv"
#connector_f <- "/scratch/trcanmed/connector/local/share/data/Models.RDs"
#connector <- read.table(connector_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
connector <- readRDS(connector_f)
connector <- connector %>% filter(Cluster %in% c(uno, due))

merged <- merge(recist, connector, by = "ShortID")

#meta_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meta_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda <- meda %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda$ShortID <- substr(meda$sample_id_R, 1, 7)
meda <- filter(meda, grepl("LMX_BASALE", type))

merged2 <- merge(merged, meda, by = "ShortID")
merged2 <- merged2[, c("sample_id_R", "ShortID", "Cluster", "batch")]
rownames(merged2) <- merged2$sample_id_R
merged2$sample_id_R <- NULL
names(merged2)[names(merged2) == "ShortID"] <- "model"

write.table(merged2, results, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)