### FOR LMX

library(tidyverse)

meda <- snakemake@input[["metadati_biobanca"]]
ne_tcf <- snakemake@input[["basali_tcf"]]
meta <- snakemake@output[["meta"]]

#ne_tcf <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/mut_samples_data_NE"
ne_tcf <- read.table(ne_tcf, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ne_tcf <- ne_tcf %>% filter(!replicates %in% c("R2", "R3"))

#meda <- "/mnt/cold1//snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meda, sep='\t', quote="", header=TRUE)
meda$sample_id <- gsub('-', '.', meda$sample_id_R, fixed = TRUE)
meda <- meda %>%
  filter(type %in% c("LMX_BASALE", "LMX_BASALE.1"))
meda$model <- substr(meda$sample_id, 1, 7)
meda <- meda[,c("model", "sample_id", "batch")]
meda <- merge(meda, ne_tcf, by="model")
meda$id <- NULL
meda$replicates <- NULL
meda$geno_mut <- gsub("NE_", "", meda$geno_mut)

# meda <- meda %>%
#   group_by(model) %>%
#   mutate(replicates = paste0("R", row_number())) %>%
#   ungroup()
meda <- as.data.frame(meda)
rownames(meda) <- meda$sample_id
meda$sample_id <- NULL

write.table(meda, meta ,quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)