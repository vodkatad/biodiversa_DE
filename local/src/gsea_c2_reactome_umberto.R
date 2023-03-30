#!/usr/bin/env Rscript

library('GSVA')
library('getopt')

train_f <- snakemake@input[["train_tmm"]]
test_f <- snakemake@input[["test_tmm"]]
rds_sign <- snakemake@input[["rds"]]
result <- snakemake@output[["gsva_res"]]

#expr_1 <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/magnum/train_tmm.tsv.gz"
#expr_2 <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/magnum/test_tmm.tsv.gz"
#rds_sign <- "/mnt/trcanmed/snaketree/prj/scRNA/dataset/CRC0327_pseudobulks/c2.symbol.rds"

train <- read.table(train_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
test <- read.table(test_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

matrix <- merge(train, test, by ="Geneid")
rownames(matrix) <- matrix$Geneid
matrix$Geneid <- NULL

geneset <- readRDS(rds_sign)
#expr_data <- read.table(gzfile(expr_file), sep="\t", header=TRUE, row.names=1)
matrix <- log(matrix+1, base=2)

#ssgsea.norm
#Barbie  et  al.   (2009)  normalizing  the  scores  by  the  absolute  difference
#between the minimum and the maximum,  as described in their paper.   Whenssgsea.norm=FALSEthis last normalization step is skipped
res <- gsva(as.matrix(matrix), geneset, kcdf="Gaussian", method="gsva")
res <- res[grepl("REACTOME", rownames(res)),]
write.table(res, file=result, quote=FALSE, sep="\t")
