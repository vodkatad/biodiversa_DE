## counts Marika singles
library(tidyverse)

c <- snakemake@input[["ck"]]
s <- snakemake@output[["samples_data"]]
wild <- snakemake@wildcards[["smodel"]]

#c <- "/mnt/cold1/bioinfotree/prj/ngs_standard_processing/dataset/RNAseq_Marika_Trametinib_03_2026/GEP.count.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(c) <- c$Geneid

cols_to_remove <- c("Geneid",
                    "CRC0322_APCmin", "CRC0322_APCplus",
                    "CRC0327_APCmin", "CRC0327_APCplus",
                    "CRC1502_NT_1", "CRC1502_TRAM_1",
                    "CRC1477_TRAM_2")

c <- c[, !(colnames(c) %in% cols_to_remove)]

id <- data.frame(id = colnames(c), stringsAsFactors = FALSE) %>%
  mutate(
    model = str_extract(id, "^[^_]+"),
    replicate = str_extract(id, "_\\d+$") %>% str_remove("_"),
    replicate = paste0("R", replicate),
    treat = str_remove(id, "^[^_]+_"),
    treat = str_remove(treat, "_\\d+$")
  ) %>%
  dplyr::select(id, model, treat, replicate)

id <- id %>% filter(model == wild)

write.table(id, file=s, row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)