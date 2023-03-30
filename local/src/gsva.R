#!/usr/bin/env Rscript

library('GSVA')
library('getopt')
opts <- matrix(c(
        'help', 'h', 0, 'logical',
        'expression' , 'e', 1, 'character',
        'signature'  , 's', 1, 'character',
        'outfile','o',1, 'character',
        'method','m',1, 'character'
        ), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$signature) || is.null(opt$expression) || is.null(opt$outfile) || is.null(opt$method)) {
    usage <- getopt(opts, usage=TRUE)
    stop(usage)
}

save.image("provagsea.R")
expr_file <- opt$expression
rds_sign <- opt$signature

geneset <- readRDS(rds_sign)
expr_data <- read.table(gzfile(expr_file), sep="\t", header=TRUE, row.names=1)
expr_data <- log(expr_data+1, base=2)

#ssgsea.norm
#Barbie  et  al.   (2009)  normalizing  the  scores  by  the  absolute  difference
#between the minimum and the maximum,  as described in their paper.   Whenssgsea.norm=FALSEthis last normalization step is skipped
res <- gsva(as.matrix(expr_data), geneset, kcdf="Gaussian", method=opt$method)
write.table(res, file=opt$outfile, quote=FALSE, sep="\t")
