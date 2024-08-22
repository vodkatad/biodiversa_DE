library(tidyverse)

metadata_o_f <- snakemake@input[["metadata"]] 
whoiswho <- snakemake@input[["buoni"]]
meta <- snakemake@output[["meta"]]

#metadata_o_f <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
tipi <- c("LMX_BASALE")
meda_f <- meda_f %>% filter(type %in% tipi)

#whoiswho <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_nolmh.tsv"
ww <- read.table(whoiswho, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(ww)[names(ww)=="CASE"] <- "sample"
res <- merge(meda_f, ww, by="sample")

## to check if the number are correct w/o replicates
#resnorep <- res
#resnorep <- resnorep[!duplicated(resnorep$sample),]
## they are ok

rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
res$type.x <- NULL
res <- res %>% filter(!type.y == "Validation not performed")

for (i in rownames(res)) {
  if (res[i, "type.y"]== "Validation successful") {
    res[i, "validation"] <- "successful"
  } else {
    res[i, "validation"] <- "notest.failed"
  }
}
res$type.y <- NULL

succ <- res %>% filter(validation == "successful")
f <- res %>% filter(validation == "notest.failed")

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)