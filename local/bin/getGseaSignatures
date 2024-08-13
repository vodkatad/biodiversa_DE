#!/usr/bin/env Rscript

library("getopt")
library("GSEABase")

opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'gmt' , 'g', 1, 'character',
        'rds'  , 'o', 1, 'character',
        'category', 'c', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$gmt) || is.null(opt$rds) || is.null(opt$category)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}

gmt_file <- opt$gmt
outputfile <- opt$rds
cat <- opt$category

genesetcollection <- getGmt(gmt_file,collectionType=BroadCollection(category=cat),geneIdType=SymbolIdentifier()) # is cat needed? FIXME
saveRDS(genesetcollection, file = outputfile)

