#!/usr/bin/env Rscript
options(error=traceback)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]

indf <- read.table(input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
indf <- indf[!is.na(indf$vol), ]
median <- median(indf$vol)
median
indf$chemio <- ifelse(indf$vol >= median, "nonresp", "resp")
table(indf$chemio)
indf$basale <- NULL
indf$type <- NULL
indf$vol <- NULL
table(indf$chemio, indf$batch)
table(indf$chemio, indf$coverage)
write.table(indf, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
