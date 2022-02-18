## prepare for add_description
## gene basemean lfc2 lfcse padj

library(tidyverse)

genes_f <- snakemake@input[["inp"]]
results <- snakemake@output[["output"]]

#genes_f <- "/scratch/trcanmed/DE_RNASeq/dataset/res_nonres_palbociclib_extreme/type_cutoff0.05-responder.vs.non_responder.deseq2.tsv"
genes_f <- read.table(genes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(genes_f)
rownames(genes_f) <- NULL
genes_f <- cbind(genes,genes_f)
genes_f <- genes_f %>% mutate(genes = gsub("H_", "", genes))
genes_f$stat <- NULL
genes_f$pvalue <- NULL

write.table(genes_f, file=results, quote = FALSE, sep = "\t", col.names = TRUE)