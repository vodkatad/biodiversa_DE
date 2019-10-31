#!/usr/bin/env Rscript
#options(error=traceback)
args <- commandArgs(trailingOnly=TRUE)
counts <- args[[1]]
metadataf <- args[[2]]
design <- args[[3]]
prefix <- args[[4]] # if all do not filter, otherwise keep ^prefix_ only genes. 
#useful for xeno mice-man counts matrices made by Ivan's pipeline.
outprefix <- args[[5]]
minc <- args[[6]]
minsamples <- args[[7]]
image <- args[[8]]
# check numeric etc

save.image(image)
library("BiocParallel")
register(MulticoreParam(4)) # TODO from CORES
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

data <- read.table(gzfile(counts), header=T, sep="\t", row=1)
if (prefix != "all") {
  data <- data[grep(paste0("^", prefix, "_"), rownames(data)),]
}

metadata <- read.table(metadataf, sep="\t", header=T, row=1)
fdesign <- as.formula(design)
print(terms(fdesign)[[2]])
rownames(metadata) <- gsub("-", ".", rownames(metadata), fixed=TRUE)
new_data <- data[,match(rownames(metadata), colnames(data))]
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=FALSE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)

save.image(image)

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

res <- data.frame(highg=highgenes, highsdgenes=highsdgenes)
write.table(res, file=paste0(outprefix, "_high.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

save.image(image)
