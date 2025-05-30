## deg DHA primo

library(tidyverse)
library(readxl)


metadata_o_f <- snakemake@input[["metadata"]] 
casi_f <- snakemake@input[["dha"]]
bb_f <- snakemake@input[["badboys"]]
meta <- snakemake@output[["sample"]]


#casi_f <- "/scratch/trcanmed/DE_RNASeq/local/share/data/DHA_primo.xlsx"
casi <- read_xlsx(casi_f)

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMO_BASALE", type))
meda_f$CASE <- substr(meda_f$sample_id_R, 1,7)

merged <- merge(casi, meda_f, by = "CASE")
res <- merged
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
names(res)[names(res) == "CASE"] <- "sample"
res$type.y <- NULL

for (i in seq(rownames(res))){
  if (res[i, "type.x"] == "responder") {
    res[i, "type"] <- "R"
  } else {
    res[i, "type"] <- "NR"
  }
}

res$type.x <- NULL

### rimozioni badboys
#bb_f <- "/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/badboys"
bb <- read.table(bb_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
bb <- bb$V1
res <- res %>% filter(!sample %in% bb)

## rimosso CRC0578
table(res$type)

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)