library(CRISclassifier, quietly=TRUE)
library("org.Hs.eg.db")

#install.packages("/home/data/work/scRNA_CRIS/ncomms15107-s20/CRISclassifier/CRISclassifier_1.0.0.tar.gz")

#crisGenesFile <- snakemake@input[["cris"]]
exprfile <- snakemake@input[["expr"]]
debug <- snakemake@params[["debug"]]
outPrefix <- snakemake@params[["prefix"]]

if (debug == "yes") {
  save.image(file=paste0(outPrefix,'.debug','.RData'))
}

#crisGenes <- read.table(crisGenesFile, sep="\t", header=TRUE)

#gbm <- load_cellranger_matrix(cellrangerDir)
#ex <- exprs(gbm)
ex <- read.table(exprfile, sep="\t", header=TRUE)

#wanted_mtx <- ex[rownames(ex) %in% ensg_cris$ENSEMBL,]
#at_least_one <- apply(wanted_mtx, 1, function(x) sum(x>0))
#n <- names(at_least_one[at_least_one>=1])
#wanted_mtx_expr <- wanted_mtx[rownames(wanted_mtx) %in% n,]
#df <- as.data.frame(as.matrix(wanted_mtx_expr))
#merged <- merge(df, ensg_cris, by.x="row.names", by.y="ENSEMBL")
#merged$ENTREZID <- NULL
#merged$Row.names <- NULL

#dims <- dim(merged)
#res <- merged[,c(dims[2],seq(1,(dims[2]-1)))]
#colnames(res)[1] <- "SYMBOL"
#tt <- table(res[,"SYMBOL"])
#remove <- names(tt[tt>1])
#res2 <- res[!res[,"SYMBOL"] %in% remove,]
ex$SYMBOL <- substr(rownames(ex), 3, nchar(rownames(ex)))
dims <- dim(ex)
res <- ex[,c(dims[2],seq(1,(dims[2]-1)))]
#colnames(res)[1] <- "SYMBOL"

write.table(res, file=paste0(outPrefix, ".tmp"), sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
cris_classifier(paste0(outPrefix, ".tmp"), outPrefix)

if (debug == "yes") {
  save.image(file=paste0(outPrefix,'.debug','.RData'))
}
