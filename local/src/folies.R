setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4')
#load('')
load('dds.Rdata')
library("BiocParallel")
metadata$typenonum <- unlist(lapply(metadata$type, function(x) {
y <- strsplit(as.character(x), '.', fixed=TRUE)[[1]][1]
# y2 <- strsplit(y, '_', fixed=TRUE)[[1]]
# last <- length(y2)
# if (!is.na(as.integer(y2[last]))) {
#   y2 <- y2[-last]
# }
# return(paste(y2, collapse="_"))
return(y)
}))
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, annotation_col=metadata[,c('typenonum','batch')],annotation_row=metadata[,c('typenonum','batch')] )
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, annotation_col=metadata[,c('typenonum','batch')],annotation_row=metadata[,c('typenonum','batch')] )
head(annotation_col=metadata[,c('typenonum','batch')])
head(metadata[,c('typenonum','batch')])
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, annotation_col=metadata[,c('typenonum','batch')] )
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, annotation_col=metadata[,c('typenonum','batch')], col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F)
dim(sampleDistMatrix)
dim(metadata)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, show_rownames=F, show_colnames=F, col=colors,annotation_row=metadata[,c('typenonum','batch')] )
metadata$typenonum
#ggplot(metadata, aes(x="", y=, fill=group))+
geom_bar(width = 1, stat = "identity")
d <- as.data.frame(table(metadata$typenonum))
ggplot(d, aes(x="", y=Freq, fill=var()))+
geom_bar(width = 1, stat = "identity")
ggplot(d, aes(x="", y=Freq, fill=Value)+
geom_bar(width = 1, stat = "identity")
)
ggplot(d, aes(x="", y=Freq, fill=Value))+
geom_bar(width = 1, stat = "identity")
head(d)
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_fill_colorblind()
#ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_fill_colorblind()
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
#ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_color_brewer(palette='Dark2')
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_color_brewer(palette='Dark2')
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_fill_brewer(palette='Dark2')
ggplot(d, aes(x="", y=Freq, fill=Var1))+
geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_bw()+scale_fill_brewer(palette='Paired')
dim(metadata)
head(fpkm)
head(fpkm_d
)
cors <- cor(fpkm_d, method="spearman")
head(cors)
pheatmap(cors, show_rownames=F, show_colnames=F,annotation_row=metadata[,c('typenonum','batch')] )
pheatmap(cors, show_rownames=F, show_colnames=F,annotation_row=metadata[,c('typenonum','batch')] )
plotPCA(vsd, intgroup='type')
plotPCA

> getMethod("plotPCA",'DESeqTransform')

ppca <- function (object, intgroup, pc1, pc2, ntop=500, returnData=FALSE) 
{

    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
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


ppca(vsd, 'type', 1, 2)
ppca(vsd, 'batch', 1, 2)
ppca(vsd, 'batch', 1, 3)
ppca(vsd, 'batch', 1, 4)
ppca(vsd, 'batch', 1, 5)
ppca(vsd, 'batch', 2, 5)
ppca(vsd, 'batch', 2, 3)

d <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/april2020/merged_hs_mm.tsv.gz', sep="\t", header=TRUE)
rownames(d) <- d$Geneid
d$Geneid <- NULL
mouse <- d[grepl("^M_",rownames(d)),]
human <- d[grepl("^H_",rownames(d)),]
plot(density(as.numeric(human)))
plot(density((human)))
plot(density(as.numeric(unlist(human))))
dim(human)
head(as.numeric(human))
head(as.numeric(as.matrix(human)))
h <- head(as.numeric(as.matrix(human)))
m <- head(as.numeric(as.matrix(m)))
m <- head(as.numeric(as.matrix(mouse)))
summary(h);
summary(m)
m <- as.numeric(as.matrix(mouse))
h <- as.numeric(as.matrix(human))
summary(m)
summary(h)
plot(density(m))
summary(h[h!=0])
plot(density(h[h!=0]))
hh <- h[h!=0]
summary(hh)
hh <- log10(h[h!=0])
mm <- log10(m[m!=0])
plot(density(hh))
plot(density(mm))
mmm <- log10(m)
hhh <- log10(h)
summary(h)
log10(0)
hhh <- log10(h+0.1)
mmm <- log10(m+0.1)
plot(density(mmm))
plot(density(hhh))
head(mouse)
mousemean <- rowMeans(mouse)
head(mousemean)
mousemean <- colMeans(mouse)
humanmean <- colMeans(human)
plot(density(mousemean))
plot(density(humanmean))
plot(density(humanmean), color="blue")+lines(density(mousemean), color="red")
plot(density(humanmean), col="blue")+lines(density(mousemean), col="red")
plot(density(mousemean), col="red")+lines(density(humanmean), col="blue")
head(mousemean)
barplot(mousemean[grepl('X',names(mousemean))],mousemean[grepl('H',names(mousemean))] )
boxplot(mousemean[grepl('X',names(mousemean))],mousemean[grepl('H',names(mousemean))] )
boxplot(mousemean[grepl('X',names(mousemean))],mousemean[grepl('LMO',names(mousemean))] )
head(mousemean)
df <- data.frame(rownames=names(mousemean), count=mousemean)
df$type <- substr(rownames(df), 0, 7)
head(df)
df$type <- substr(rownames(df), 7,3)
head(df)
df$type <- substr(rownames(df), 7,10)
head(df)
df$type <- substr(rownames(df), 8,10)
head(df)
ggplot(ToothGrowth, aes(x=dose, y=len)) +
geom_boxplot()
ggplot(d, aes(x=type, y=count)) + geom_boxplot()
ggplot(df, aes(x=type, y=count)) + geom_boxplot()
ggplot(df, aes(x=type, y=count)) + geom_violin()
ggplot(df, aes(x=type, y=count)) + geom_violin()
ggplot(df, aes(x=type, y=count)) + geom_boxplot()
ggplot(df, aes(x=type, y=count, fill=count)) + geom_boxplot()
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()
library(ggsignif)
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE)
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMO','LMX')))
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMX')))
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMX'), c("LMO","LMX"),c('PRH','PRX')))
head(d)
head(df)
df$hx <- ifelse(df$type %in% c("LMX",'PMX','PRX'), 'X', 'H')
ggplot(df, aes(x=hx, y=count, fill=hx)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("X","H"))
)
wilcox.test(df[df$type=="X",'count'],df[df$type=="H",'count'])
wilcox.test(df[df$hx=="X",'count'],df[df$hx=="H",'count'])
t.test(df[df$hx=="X",'count'],df[df$hx=="H",'count'])

setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4_cetuxi_treat_PDX')
load('dds.Rdata')
library(DESeq2)
cetuxi <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/april2020/selected_metadata', header=F, sep="\t")
p <- plotPCA(vsd, 'treat', returnData=T)
m <- merge(p, cetuxi, by.x="row.names", by.y="V1")
m$recist <- ifelse(m$V4 < -0.50, 'OR', ifelse(m$V4 > 0.35, 'PD', 'SD'))
ggplot(data=m, aes(PC1, PC2, color=as.factor(recist)))+geom_point()+theme_bw()

