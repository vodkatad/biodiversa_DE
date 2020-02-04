#!/usr/bin/env Rscript
#options(error=traceback)
args <- commandArgs(trailingOnly=TRUE)
counts <- args[[1]]
len <- args[[2]]
prefix <- args[[3]] # if all do not filter, otherwise keep ^prefix_ only genes. 
#useful for xeno mice-man counts matrices made by Ivan's pipeline.
minc <- as.numeric(args[[4]])
minsamples <- as.numeric(args[[5]])
image <- args[[6]]
fpkmf <- args[[7]]
outprefix <- args[[8]]
design <- "~1"
fdesign <- as.formula(design)
save.image('pippo.Rdata')

library("BiocParallel")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")

register(MulticoreParam(4)) # TODO from CORES
library(DESeq2)

data <- read.table(gzfile(counts), header=T, sep="\t", row=1)
if (prefix != "all") {
  data <- data[grep(paste0("^", prefix, "_"), rownames(data)),]
}

dds <- DESeqDataSetFromMatrix(countData = data, colData=data.frame(row.names=colnames(data)), design=fdesign)
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=FALSE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatf <- paste0(outprefix, "_vsd_dists.pdf")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, file=pheatf)
if (design != "~1") {
    out <- tryCatch(
        {
           ints <- attr(terms(fdesign),"term.labels")
           for (i in seq(1, length(ints))) {
               plotPCA(vsd, intgroup=ints[i])
               ggsave(paste0(outprefix, "_", ints[i], "_pca.pdf"))
           }
        },
        error=function(cond) {
            message(cond)
        },
        warning=function(cond) {
            message(cond)
        },
        finally=function(cond) { dev.off() }
    )
}

cc <- counts(dds,normalized=TRUE)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
if (design != "~1") {
    df <- as.data.frame(colData(dds)[,attr(terms(fdesign),"term.labels"), drop=F])
} else {
    df <- as.data.frame(colData(dds)[,c(terms(fdesign)[[2]]),drop=F])
}
pheatf2 <- paste0(outprefix, "_cc_high.pdf")
pheatmap(cc[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, file=pheatf2)
highgenes <- rownames(cc[select,])
sds <- apply(cc, 1, sd)
highsd <- cc[order(sds, decreasing=TRUE)[1:20],]
pheatf3 <- paste0(outprefix, "_cc_highsd.pdf")
pheatmap(highsd, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, file=pheatf3)
highsdgenes <- rownames(highsd)

lens <- read.table(len, sep="\t", header=T)
order <- rownames(mcols(dds))
lens <- lens[match(order, lens$Geneid), ]
all(lens$Geneid == order)
mcols(dds)$basepairs  <- lens$length
fpkm_d <- fpkm(dds)
write.table(fpkm_d, gzfile(fpkmf), quote=F, row.names=T, col.names=T, sep="\t")

save.image(image)
