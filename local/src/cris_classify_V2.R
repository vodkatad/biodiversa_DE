library(CRISclassifier, quietly=TRUE)


exprfile <- snakemake@input[["expr"]]
outPrefix <- snakemake@params[["prefix"]]

ex <- read.table(exprfile, sep="\t", header=TRUE)

ex$SYMBOL <- rownames(ex)
dims <- dim(ex)
res <- ex[,c(dims[2],seq(1,(dims[2]-1)))]

write.table(res, file=paste0(outPrefix, ".tmp"), sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
cris_classifier(paste0(outPrefix, ".tmp"), outPrefix)
