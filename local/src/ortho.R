### TODO maybe more generic for the future?

#!/usr/bin/env Rscript

library(biomaRt)

set.seed(42)

data <- snakemake@input[[1]]
output <- snakemake@output[[1]]


human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genes <- read.table(data, sep='\t', quote="", header=TRUE)


ortho <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                values = genes$Gene, mart=mouse,
                attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

df <- merge(ortho, genes, by.x='HGNC.symbol', by.y='Gene')
df[,1] <- NULL
names(df)[1] <- c("Gene")

write.table(df, output, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
