
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
meta$sample_id_R <- gsub('-','.', meta$sample_id_R, fixed=TRUE)
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

### Francesco's markers
pdata2 <- data.frame(HPRT1=as.numeric(data[rownames(data)=="H_HPRT1",]), CETN2=as.numeric(data[rownames(data)=="H_CETN2",]))
mpdata2 <- melt(pdata2)
ggplot(data=mpdata2, aes(x=value, fill=variable))+geom_histogram(position="dodge")+theme_bw()+facet_wrap(~variable)

meta <- meta[match(colnames(data), meta$sample_id_R), ]
pdata2 <- data.frame(HPRT1=as.numeric(data[rownames(data)=="H_HPRT1",]), CETN2=as.numeric(data[rownames(data)=="H_CETN2",]),CTCF=as.numeric(data[rownames(data)=="H_CTCF",]), ATOH1=as.numeric(data[rownames(data)=="H_ATOH1",]), type=meta$type)
pdata2$sample <- as.factor(sapply(strsplit(as.character(pdata2$type), ".", fixed=TRUE), function(x) {x[1]}))
pdata2 <- pdata2[order(pdata2$sample),]


mpdata2 <- melt(pdata2)


single_plot <- function(data, gene) {
 d <- data[data$variable==gene,]
 d$x = seq(1, nrow(d))
 print(ggplot(data=d, aes(y=value,x=sample,fill=sample))+geom_jitter(position=position_jitter(0.2))+geom_violin(alpha=0.4,trim=FALSE)+geom_boxplot(width=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 15))+ggtitle(gene))

###

}



single_plot_a <- function(data, gene) {
d <- data[data$variable==gene,]
d$x = seq(1, nrow(d))
print(ggplot(data=d, aes(y=value,x=x,fill=sample))+geom_bar(stat="identity")+theme_bw()+ggtitle(gene))
}

get_sd_deciles <- function(data) {
  s <- seq(0,1, length=101)
  data$deciles <- cut(data$sd, quantile(data$sd, prob = s, type = 5), include.lowest=TRUE)
  levels(data$deciles) <- s
  return(data)
}


single_plot_enh <- function(data, meta, gene, log=FALSE) {
  pdata2 <- data.frame(gene=as.numeric(data[rownames(data)==gene,]), type=meta$type)
  pdata2$sample <- as.factor(sapply(strsplit(as.character(pdata2$type), ".", fixed=TRUE), function(x) {x[1]}))
  pdata2 <- pdata2[order(pdata2$sample),]
  d <- melt(pdata2)
  if (log) {
    d$value <- log(d$value)/log(2)
  }
  
  d$x = seq(1, nrow(d))
  print(ggplot(data=d, aes(y=value,x=sample,fill=sample))+geom_jitter(position=position_jitter(0.2))+geom_violin(alpha=0.4,trim=FALSE)+geom_boxplot(width=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 15))+ggtitle(gene))
  
  ###
  
}

### FPKM
data <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/fpkm.tsv.gz'), sep="\t", header=TRUE, row.names=1)

meta <- read.table('/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020/selected_metadata_annot_final_nolinfo_nooutlier', sep="\t", header=TRUE)
meta$sample_id_R <- gsub('-','.', meta$sample_id_R, fixed=TRUE)
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

high <- getlowestdecilesd(pdata[pdata$meandeciles==0.9,])
gethighestdecilesd(pdata[pdata$meandeciles==0.9,])

single_plot_enh_s <- function(data, meta, gene, log=FALSE) {
  pdata2 <- data.frame(gene=as.numeric(data[rownames(data)==gene,]), type=meta$type)
  pdata2$sample <- as.factor(sapply(strsplit(as.character(pdata2$type), ".", fixed=TRUE), function(x) {x[1]}))
  pdata2 <- pdata2[order(pdata2$sample),]
  d <- melt(pdata2)
  if (log) {
    d$value <- log(d$value)/log(2)
  }
  
  d$x = seq(1, nrow(d))
  ggplot(data=d, aes(y=value,x=sample,fill=sample))+geom_jitter(position=position_jitter(0.2))+geom_violin(alpha=0.4,trim=FALSE)+geom_boxplot(width=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 15))+ggtitle(gene)
  ggsave(paste0(gene ,'violin.png'))
  ###
}

############


###############
library(ggplot2)
starok <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK/vsd.tsv.gz', sep="\t", header=T)
starko <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/vsd.tsv.gz', sep="\t", header=T)
g_sok <- rownames(starok)
g_sko <- rownames(starko)
g <- intersect(g_sok, g_sko)
d_starok <- starok[rownames(starok) %in% g,]
d_starko <- starko[rownames(starko) %in% g,]
dim(starko)
dim(starok)
length(g)
length(g_sok)
length(g_sko)
dim(d_starok)
dim(d_starko)
all(rownames(d_starko)==rownames(d_starko))
all(colnames(d_starko)==colnames(d_starko))
correlations_samples <- cor(d_starko, d_starok)
correlations_genes <- cor(t(d_starko), t(d_starok))

bp <- function(cors, title) {
  matched <- diag(cors)
  other <- cors[lower.tri(cors)]
  pd <- data.frame(values=c(matched, other), class=c(rep('diagonal',length(matched)), rep('other', length(other))))
  print(ggplot(data=pd, aes(y=values, x=class, color=class))+geom_boxplot()+theme_bw()+ggtitle(title)+theme(text=element_text(size=21)))
}

cor_genes <- diag(correlations_genes)
other_genes <- correlation_genes[lower.tri(correlations_genes)]

oth <- other_genes[sample.int(581285656, size =34097)]
pd <- data.frame(values=c(cor_genes, oth), class=c(rep('diagonal',length(cor_genes)), rep('other', length(oth))))
ggplot(data=pd, aes(y=values, x=class, color=class))+geom_boxplot()+theme_bw()+ggtitle('genes')+theme(text=element_text(size=21))



