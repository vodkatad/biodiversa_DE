## counts KRAS
library(tidyverse)

c <- snakemake@input[["ck"]]
s <- snakemake@output[["samples_data"]]
wild <- snakemake@wildcards[["smodel"]]

#c <- "/mnt/cold1/bioinfotree/prj/ngs_standard_processing/dataset/KRASinh_simo/GEP.count.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

id <- colnames(c)[2:25]
id <- as.data.frame(id)
id[10, "id"] <- "CRC0031_M_C_T_EXP2"

id <- id %>%
  mutate(
    model = str_extract(id, "^[^_]+"),
    replicate = str_extract(id, "EXP\\d+$"),
    replicate = ifelse(is.na(replicate), "EXP1", replicate),
    treat = str_remove(id, "^[^_]+_"),
    treat = str_remove(treat, "_EXP\\d+$")
  ) %>%
  select(id, model, treat, replicate)

id[10, "id"] <- "CRC0031_M_C_TEXP2"

id <- id %>% filter(model == wild)

write.table(id, file=s, row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)
