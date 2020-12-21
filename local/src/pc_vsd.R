
ppca <- function (data, meta, intgroup, pc1, pc2, ntop=500, returnData=FALSE)
{
  
  rv <- rowVars(data)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(data[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(meta))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(meta[, intgroup, 
                                               drop = FALSE])
  group <- meta[, intgroup]
  d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group,
                 name = colnames(data))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  d$group <- as.factor(d$group) 
  print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 2) + xlab(paste0("PC", pc1, ': ', round(percentVar[1] *  
                                                                100), "% variance")) + ylab(paste0("PC", pc2, ":", round(percentVar[2] * 
                                                                                                                           100), "% variance"))+theme_bw())
  return(percentVar)
}


cppca <- function (data, meta, intgroup, pc1, pc2, ntop=500, returnData=FALSE)
{
  
  rv <- rowVars(data)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(data[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(meta))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(meta[, intgroup, 
                                    drop = FALSE])
  group <- meta[, intgroup]
  d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group,
                  name = colnames(data))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
          geom_point(size = 2) + xlab(paste0("PC", pc1, ': ', round(percentVar[1] *  
                                                                      100), "% variance")) + ylab(paste0("PC", pc2, ":", round(percentVar[2] * 
                                                                                                                                 100), "% variance")) +theme_bw())
  return(pca$rotation)
}

metadata <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/magnum/selected_metadata.tsv', sep="\t", header=TRUE, row.names = 1)
data <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/magnum/selected_matrix.tsv.gz'), sep="\t", header=TRUE, row.names=1)
metadata <- metadata[match(colnames(data), rownames(metadata)),]
all(rownames(metadata)== colnames(data))
ppca(as.matrix(data), metadata, 'batch', 1,2)
ppca(as.matrix(data), metadata, 'batch', 2,3)
ppca(as.matrix(data), metadata, 'batch', 1,3)
ppca(as.matrix(data), metadata, 'batch', 1,2)
# PC4 (7% variance) sees something between batch 4-5 / 1-2

library(limma)
mat <- limma::removeBatchEffect(data, metadata$batch)
ppca(as.matrix(mat), metadata, 'batch', 2,4)
ppca(as.matrix(mat), metadata, 'batch', 1,2)

##???

###################
# TODO REMOVE lymphomas and outlier!
data <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/vsd.tsv.gz'), sep="\t", header=TRUE, row.names=1)

meta <- read.table('/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020/selected_metadata_annot_final_nolinfo_nooutlier', sep="\t", header=TRUE)
data <- data[, colnames(data) %in% meta$sample_id_R,]

sds <- apply(data, 1, sd)
means <- rowMeans(data)
library(ggplot2)
pdata <- data.frame(row.names=names(means), mean=means, sd=sds, median=apply(data, 1, median))
library(reshape2)
mpdata <- melt(pdata)
#library(DESeq2)
#load('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/dds.Rdata')
#library(vsn)
#meanSdPlot(assay(vsd))
#pdata$mrank <- rank(pdata$mean)
#plot(pdata$mrank, pdata$sd)


# for each decile of expression we extract the lowest decile in sd inside that selection of genes to obtain housekeeping genes
pdata$meandeciles <- cut( pdata$mean, quantile(pdata$mean, prob = seq(0, 1, length = 11), type = 5), include.lowest=TRUE )
levels(pdata$meandeciles) <- seq(0,1, length=11)
#pdata$md <- as.character(pdata$meandeciles)
getlowestdecilesd <- function(data) {
  s <- seq(0,1, length=101)
  data$deciles <- cut(data$sd, quantile(data$sd, prob = s, type = 5), include.lowest=TRUE)
  levels(data$deciles) <- s
  return(rownames(data[data$deciles == s[1],]))
}

gethighestdecilesd <- function(data) {
  s <- seq(0,1, length=101)
  data$deciles <- cut(data$sd, quantile(data$sd, prob = s, type = 5), include.lowest=TRUE)
  levels(data$deciles) <- s
  return(rownames(data[data$deciles == s[length(s)-1],]))
}

getlowestdecilesd(pdata[pdata$meandeciles==0.9,])
gethighestdecilesd(pdata[pdata$meandeciles==0.9,])

#hist(as.numeric(data[rownames(data)=="H_CTCF",]))
pdata2 <- data.frame(CTCF=as.numeric(data[rownames(data)=="H_CTCF",]), IGFBP2=as.numeric(data[rownames(data)=="H_IGFBP2",]), LGALS3=as.numeric(data[rownames(data)=="H_LGALS3",]))
mpdata2 <- melt(pdata2)
ggplot(data=mpdata2, aes(x=value, fill=variable))+geom_histogram(position="dodge")+theme_bw()


s <- seq(0,1, length=101)
pdata$deciles <- cut(pdata$sd, quantile(pdata$sd, prob = s, type = 5), include.lowest=TRUE)
levels(pdata$deciles) <- s
pdata[rownames(pdata)=="H_LGALS3",]

