#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
matrix_f <- args[[1]]
gct_f <- args[[2]]

d <- read.table(matrix_f, sep="\t", stringsAsFactors=FALSE, row.names=NULL)

d$Description <- rep("NA", nrow(d))
d <- d[, c(1, ncol(d), seq(3, ncol(d)-1))]
colnames(d)[1] <- 'NAME'

write.table(d, file=gct_f, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE)
