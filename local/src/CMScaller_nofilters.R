### CMScaller

library(tidyverse)
library(CMScaller)

expr <- snakemake@input[["expr"]]
plot <- snakemake@output[["CMS_heatmap"]]
results <- snakemake@output[["RES"]]

expr <- read.table(expr, quote = "", sep = "\t", header = TRUE)
expr <- expr %>% remove_rownames %>% column_to_rownames(var="entrez")

pdf(plot)
res <- CMScaller(emat=expr, RNAseq=TRUE, FDR=0.05)
dev.off()

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
