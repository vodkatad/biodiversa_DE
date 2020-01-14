#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]


#CRC0018LMX0B02204TUMR14000      CRC0018LMX      -0.56   first   LMX     low     basale

vols <- read.table(input, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(vols) <- c("long_gid", "gid", "perc","batch","type","coverage","basale")

vols$class <- ifelse(vols$perc < -0.5, 'OR', ifelse(vols$perc > 0.35, 'PD', 'SD'))

table(vols$class)
table(vols$class, vols$batch)
table(vols$class, vols$coverage)


write.table(vols, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
