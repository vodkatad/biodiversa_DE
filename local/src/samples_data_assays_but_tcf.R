library(tidyverse)

sbcat <- snakemake@input[["assay"]]
s <- snakemake@input[["samples"]]
r <- snakemake@output[["meta"]]

sbcat <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/b2_general/samples_data"
sbcat <- read.table(sbcat, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sbcat <- sbcat %>% filter(geno == "b2")
sbcat <- unique(sbcat$model)
sbcat <- gsub("sh_", "", sbcat)

s <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/general/table_cases.tsv"
s <- read.table(s, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
s <- s %>% filter(model %in% sbcat)

res <- s

write.table(res, file=r, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)