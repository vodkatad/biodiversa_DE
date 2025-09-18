library(tidyverse)

ss <- snakemake@input[["meta"]]
c <- snakemake@input[["counts"]]
cf <- snakemake@output[["fcounts"]]
ssm <- snakemake@output[["meta_filtered"]]

#ss <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCGA/samples_data"
samples_data <- read.table(ss, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#c <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/ilio_wgd/TCGA_COAD-READ_counts_v3.tsv.gz"
counts <- read.table(c, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

sample_ids <- samples_data$Sample.ID
filtered_counts <- counts[, colnames(counts) %in% sample_ids, drop = FALSE]

samples_data <- samples_data %>% filter(Sample.ID %in% colnames(filtered_counts))

write.table(filtered_counts, cf, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(samples_data, ssm, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

c <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/ilio_wgd/TCGA_COAD-READ_counts.tsv.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

c <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/merged_hs_mm.tsv.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
