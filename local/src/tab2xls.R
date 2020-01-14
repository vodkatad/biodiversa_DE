#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]

library(WriteXLS)
d <- read.table(input, sep="\t", header=T, comment.char="")
WriteXLS(d, ExcelFileName=output, row.names=FALSE, col.names=TRUE)
