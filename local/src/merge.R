#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
output <- args[[1]]

merged <- read.table(gzfile(args[[2]]), sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
for (i in 3:length(args)) {
    nextdf <- read.table(gzfile(args[[i]]), sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
    merged <- merge(merged, nextdf, by="row.names")
    rownames(merged) <- merged$Row.names
    merged$Row.names <- NULL
}
merged$Geneid <- rownames(merged)
n <- ncol(merged)
merged <- merged[, c(n,seq(1, n-1))]
write.table(merged, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
