## counts totali plus counts primo

library(tidyverse)

counts_primo <- snakemake@input[["primo"]]
counts_totali <- snakemake@input[["totali"]]
counts_finali <- snakemake@output[["counts"]]

#counts_primo <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/scTimeCourse_Palbo/counts.tsv.gz"
cp <- read.table(gzfile(counts_primo), quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cp <- cp[c("CRC1257", "CRC1448")]
cp$Geneid <- rownames(cp)
rownames(cp) <- NULL
for (i in seq(length(cp$Geneid))) {
  cp[i,"Geneid"] <- paste0("H_", cp[i,"Geneid"])
}

#counts_totali <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/merged_hs_mm.tsv.gz"
ct <- read.table(gzfile(counts_totali), quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts <- ct
counts$type <- substr(counts$Geneid, 1, 1)
counts <- counts %>% filter(!type == "M")
counts$type <- NULL

merged_counts <- merge(counts, cp, by = "Geneid")
rownames(merged_counts) <- merged_counts$Geneid
merged_counts$Geneid <- NULL

write.table(merged_counts, counts_finali, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
