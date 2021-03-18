#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)

set.seed(42)

data <- snakemake@input[[1]]
output <- snakemake@output[[1]]

# setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/')
load(data)
ppca <- function (object, intgroup, pc1, pc2, ntop=500, returnData=FALSE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
    length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group,
    intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    d$group <- as.factor(d$group)
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 2) + xlab(paste0("PC", pc1, ': ', round(percentVar[1] *
    100), "% variance")) + ylab(paste0("PC", pc2, ":", round(percentVar[2] *
    100), "% variance")) + coord_fixed()+theme_bw()
}

pcs <- ppca(vsd, 'type', 1, 2, returnData=TRUE)
rownames(pcs) <- gsub('.', '-', rownames(pcs), fixed = TRUE) # solito problema del . nei replicati che in altri file è un -


write.table(pcs, file=output, sep="\t", quote=F)
