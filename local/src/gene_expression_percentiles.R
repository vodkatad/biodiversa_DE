data <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5/vsd.tsv.gz'), sep="\t", header=TRUE, row.names=1)
#/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK/
data <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz'), sep="\t", header=TRUE, row.names=1)

meta <- read.table('/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020/selected_metadata_annot_final_nolinfo_nooutlier_replisafe', sep="\t", header=TRUE)
#meta$sample_id_R <- gsub('-','.', meta$sample_id_R, fixed=TRUE)
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
