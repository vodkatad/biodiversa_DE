## correct counts ilio 

library(tidyverse)
library(org.Hs.eg.db)

c <- snakemake@input[["counts_ilio"]]
cc <- snakemake@output[["filtered_counts"]]

#c <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/ilio_wgd/TCGA_COAD-READ_counts.tsv.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ids_clean <- c$Ensembl_ID
ids_clean <- sub("\\..*", "", ids_clean)

symbols <- mapIds(org.Hs.eg.db,
                  keys = ids_clean,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
c$Gene <- symbols

analysis <- c
analysis <- analysis[,c("Ensembl_ID", "Gene")]
analysis$checknan <- grepl("nan", analysis$Gene)
analysis <- analysis %>%
  group_by(Gene) %>%
  mutate(checkdupli = n() > 1) %>%
  ungroup()

analysis <- analysis %>% filter(!checkdupli == TRUE)
analysis$checknan <- NULL
analysis$checkdupli <- NULL

counts <- merge(c, analysis, by="Ensembl_ID")
counts$Ensembl_ID <- NULL
rownames(counts) <- counts$Gene.x
counts$Gene.x <- NULL
counts$Gene.y <- NULL

#table(apply(counts, 2, is.integer))
forceMatrixToInteger <- function(m){
  apply (m, c (1, 2), function (x) {
    (as.integer(x))
  })
}

counts <- forceMatrixToInteger(counts)

write.table(counts, cc, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)