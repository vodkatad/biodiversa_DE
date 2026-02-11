## counts KRAS
library(tidyverse)

c <- snakemake@input[["counts"]]
ck <- snakemake@output[["counts_kras"]]

#c <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/counts_kras.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

c <- c %>%
  select(matches("Geneid|CRC0031|CRC1598"))

colnames(c) <- gsub("\\.", "_", colnames(c))

write.table(c, file=ck, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)