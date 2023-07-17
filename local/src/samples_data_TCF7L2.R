library(tidyverse)
library(readxl)

samples_f <- snakemake@input[["samples_or"]] 
meta <- snakemake@output[["meta"]]
case <- snakemake@wildcards[['smodel']]

#samples <- read_xlsx("/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/TCF7L2_2nd/samples_names.xlsx")
samples <- read_xlsx(samples_f)
samples$...4 <- NULL
samples$...5 <- NULL
samples <- samples[-c(133:171),]
samples$id <- samples$REPLICATES
samples$model <- substr(samples$REPLICATES, 1,7)
samples$geno <- substr(samples$REPLICATES, 9, 10)
samples$replicates <- substr(samples$REPLICATES, 12, 13)
samples$sample <- NULL
samples$`tube name` <- NULL
samples$REPLICATES <- NULL

#samples <- samples %>% filter(model == case)
write.table(samples, file=meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)