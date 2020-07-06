dd$type <- unlist(lapply(dd$type, function(x) {
y <- strsplit(as.character(x), '.', fixed=TRUE)[[1]][1]
# y2 <- strsplit(y, '_', fixed=TRUE)[[1]]
# last <- length(y2)
# if (!is.na(as.integer(y2[last]))) {
#   y2 <- y2[-last]
# }
# return(paste(y2, collapse="_"))
return(y)
}))
geom_boxplot() +theme_bw()
ggplot(ddd, aes(type, DEFA5, fill = type)) +
############################# whole or partial samples for DE
counts <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/april2020/merged_hs_mm.tsv.gz"
data <- read.table(gzfile(counts), header=T, sep="\t", row=1)
if (prefix != "all") {
data <- data[grep(paste0("^", prefix, "_"), rownames(data)),]
}
head(data)
new_data <- data[,match(rownames(metadata), colnames(data))]
head(rownames(metadta))
head(rownames(metadata))
head(colnames(data))
dim(data)
dim(metadata)
length(intersect(rownames(metadata), colnames(data)))
setdiff(rownames(metadata), colnames(data))
rownames(metadata) <- gsub("-", ".", rownames(metadata), fixed=TRUE)
new_data <- data[,match(rownames(metadata), colnames(data))]
setdiff(rownames(metadata), colnames(data))
setdiff(colnames(data), rownames(metadata))
colnames(data) <- gsub("-2", "", colnames(data), fixed=TRUE)
setdiff(colnames(data), rownames(metadata))
colnames(data) <- gsub(".2", "", colnames(data), fixed=TRUE)
setdiff(colnames(data), rownames(metadata))
new_data <- data[,match(rownames(metadata), colnames(data))]
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
head(metadata()
head(metadata)
head(metadata)
metadata$batch <- as.factor(metadata$batch)
metadata$time <- as.factor(metadata$time)
head(metadata)
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
res <- result(dds, alpha=0.05, contrast=c("treat","cetuxi","NT"), parallel=TRUE)
res <- results(dds, alpha=0.05, contrast=c("treat","cetuxi","NT"), parallel=TRUE)
head(res)
plot_volcano <- function(resnona, alpha, lfc, outfile, title) {
#title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
#resnona <- res[!is.na(res$padj),]
resnona$sign <- ifelse(abs(resnona$log2FoldChange) > lfc & resnona$padj < alpha, "both", ifelse(abs(resnona$log2FoldChange) > lfc, "LFC",
ifelse(resnona$padj < alpha, "padj", "NS")))
resnona$sign <- factor(resnona$sig, levels=c("LFC", "padj", "both", "NS"))
resnona[resnona$padj ==0,"padj"] <- .Machine$double.xmin
p <- ggplot(resnona, aes(log2FoldChange, -log10(padj))) +
geom_point(aes(col = sign),size=0.5) + theme_bw() +
scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE) + # red orange green black -> orange blue green gray
ggtitle(title)
nsign <- nrow(resnona[resnona$sig=="both",])
if (nsign > 20) {
p + geom_text_repel(data=resnona[1:10,], aes(label=rownames(resnona)[1:10]))
} else {
p + geom_text_repel(data=resnona[resnona$sig=="both",], aes(label=rownames(resnona[resnona$sig=="both",])))
}
}
resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])
title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
tile
title
plot_volcano(resnona_df, alpha, lfc, volcano, title)
plot_volcano(resnona_df, alpha, 0.5849625, volcano, title)
plot_volcano(resnona_df, 0.05, 0.5849625, volcano, title)
library(ggrepel)
plot_volcano(resnona_df, 0.05, 0.5849625, volcano, title)
write.table(resnona_df, file="test_all_cetuxi.tsv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
res2 <- read.table("../Biodiversa_up4_cetuxi_treat_PDX/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv", sep="\t", header=TRUE)
head(Res2)
head(res2)
mm <- merge(res2, res, by="row.names")
head(res)
mm <- merge(res2, resnona_df, by="row.names")
dim(mm)
dim(res2)
dim(resnona_df)
head(mm)
plot(mm$log2FoldChange.x, mm$log2FoldChange.x)
plot(mm$log2FoldChange.x, mm$log2FoldChange.y)
cor.test(mm$log2FoldChange.x, mm$log2FoldChange.y)
cor.test(mm$pvalue.y, mm$pvalue.x)
cor.test(mm$pvalue.y, mm$pvalue.y)
plot(mm$pvalue.x, mm$pvalue.y)
plot(log10(mm$pvalue.x), log10(mm$pvalue.y))
cor.test(log10(mm$pvalue.y), log10(mm$pvalue.y))
cor.test(log10(mm$pvalue.x), log10(mm$pvalue.y))
head(res)
head(res2)
sign_partial <- res[res$padj<0.05 & abs(res$log2FoldChange) > 0.58,]
sign_partial <- res2[res2$padj<0.05 & abs(res2$log2FoldChange) > 0.58,]
sign_whole <- res2[resnona_df$padj<0.05 & abs(resnona_df$log2FoldChange) > 0.58,]
dim(sign_partial)
dim(sign_whole)
length(intersect(rownames(sign_partial), rownames(sign_whole)))
library(VennDiagram)
plot(mm$lfcSE.x, mm$lfcSE.y)
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
set1 <- rownames(sign_partial)
set2 <- rownames(sign_whole)
venn.diagram(
x = list(set1, set2),
category.names = c("Set 1" , "Set 2 " , "Set 3"),
output=TRUE,
# Output features
imagetype="png" ,
height = 480 ,
width = 480 ,
resolution = 300,
compression = "lzw",
# Circles
lwd = 2,
lty = 'blank',
fill = myCol,
# Numbers
cex = .6,
fontface = "bold",
fontfamily = "sans",
# Set names
cat.cex = 0.6,
cat.fontface = "bold",
cat.default.pos = "outer",
cat.pos = c(-27, 27, 135),
cat.dist = c(0.055, 0.055, 0.085),
cat.fontfamily = "sans",
rotation = 1
)
venn.diagram(
x = list(set1, set2),
category.names = c("Set 1" , "Set 2 " , "Set 3"),
output=TRUE,
# Output features
imagetype="png" ,
height = 480 ,
width = 480 ,
resolution = 300,
compression = "lzw",
# Circles
lwd = 2,
lty = 'blank',
fill = myCol,
# Numbers
cex = .6,
fontface = "bold",
fontfamily = "sans",
# Set names
cat.cex = 0.6,
cat.fontface = "bold",
cat.default.pos = "outer",
cat.pos = c(-27, 27, 135),
cat.dist = c(0.055, 0.055, 0.085),
cat.fontfamily = "sans",
rotation = 1, filename="venn_lfc.png"
)
venn.diagram(
x = list(set1, set2),
category.names = c("Set 1" , "Set 2 " , "Set 3"),
output=TRUE,
# Output features
imagetype="png" ,
height = 480 ,
width = 480 ,
resolution = 300,
compression = "lzw",
filename="venn_lfc.png"
)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
height = 480 ,
width = 480 ,
resolution = 300,
compression = "lzw",
filename="venn_lfc.png"
)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png"
)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
sign_whole <- res2[resnona_df$padj<0.05,]
sign_partial <- res[res$padj<0.05 ,]
sign_partial <- resnona_df[resnona_df$padj<0.05 ,]
set1 <- rownames(sign_partial)
set2 <- rownames(sign_whole)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
sign_partial <- resnona_df[resnona_df$padj<0.0001 ,]
sign_whole <- resnona_df[resnona_df$padj<0.0001 ,]
set1 <- rownames(sign_partial)
set2 <- rownames(sign_whole)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
getwd()
length(set1)
sign_whole <- res2[resnona_df$padj<0.001 & abs(resnona_df$log2FoldChange) > 1,]
sign_partial <- res2[res2$padj<0.001 & abs(res2$log2FoldChange) > 1,]
set1 <- rownames(sign_partial)
set2 <- rownames(sign_whole)
venn.diagram(
x = list(set1, set2),
category.names = c("Partial" , "Whole"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
###############33 fpkm and defa5
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4/')
library(DESeq2)
load('dds.Rdata')
d <- fpkm_d
plot(density(d))
min(d[d!=0])
ld <- log(d+0.00001)
plot(density(ld))
#defa5 <- data.frame(value=ld)
head(ld)
defa5 <- data.frame(value=ld[rownames(ld)=="H_DEFA5", ], type=metadata$type)
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
defa5 <- data.frame(value=ld[rownames(ld)=="H_DEFA5", ], type=metadata$typenonum)
library(ggplot)
library(ggplot2)
ggplot(defa5, aes(x=type, y=value)) + geom_violin()+theme_bw()
ggplot(defa5, aes(x=type, y=value, fill=type())) + geom_violin()+theme_bw()
ggplot(defa5, aes(x=type, y=value, fill=type) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
)
ggplot(defa5, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun.y=median, geom="point", size=2, color="red")
ggplot(defa5, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun.y=median, geom="point", shape=23, size=2)
pd <- data.frame(value=ld[rownames(ld)=="H_TP53INP2", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)
pd <- data.frame(value=ld[rownames(ld)=="H_SOX4", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)
pd <- data.frame(value=ld[rownames(ld)=="H_FGFBP1", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)
pd <- data.frame(value=ld[rownames(ld)=="H_TOP1", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)+ggtitle('TOP1')
pd <- data.frame(value=ld[rownames(ld)=="H_SOX4", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)+ggtitle('SOX4')
pd <- data.frame(value=ld[rownames(ld)=="H_DNMBP", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)+ggtitle('DNMBP')
plot(density(ld))
plot(density(ld[ld > -10.5]))
dim(d[d==0]()
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
plot(density(apply(ld, 1, mean)))
pd <- data.frame(value=ld[rownames(ld)=="H_CMET", ], type=metadata$typenonum)
pd <- data.frame(value=ld[rownames(ld)=="H_MET", ], type=metadata$typenonum)
ggplot(pd, aes(x=type, y=value, fill=type)) + geom_violin()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ stat_summary(fun=median, geom="point", shape=23, size=2)+ggtitle('DNMBP')
pd <- data.frame(value=ld[rownames(ld)=="H_DEFA5", ], type=metadata$typenonum)
plot(density(pd))
plot(density(pd$value))
histogram(pd$value))
histogram(pd$value)
ggplot(data=pd, aes(value)) +  geom_histogram()+theme_bw()
ggplot(data=pd, aes(value), fill="darkgreen") +  geom_histogram()+theme_bw()
ggplot(data=pd, aes(value, fill="darkgreen")) +  geom_histogram()+theme_bw()
ggplot(data=pd, aes(value)) +  geom_histogram(fill="darkgreen")+theme_bw()
d <- read.table('~/matrice_microarray.tsv.gz', sep="\t", header=T)
head(d)
pd <- data.frame(value=ld[rownames(ld)=="DEFA5", ])
pd <- data.frame(value=d[rownames(d)=="DEFA5", ])
plot(density(d))
head(d)
dim(d)
plot(density(as.numeric(d)))
class(d)
head(colnames(d))
head(rownames(d))
class(d$CRC0014LM)
class(d$CRC0018LM.1)
apply(d, 1, class)
nn <- apply(d, 1, class)
nn[nn!="numeric"]
plot(density(d))
plot(density(apply(d,1, mean))
)
dim(d)
#ggplot(data=pd, aes(value)) +  geom_histogram(fill="darkgreen")+theme_bw()
pd <- data.frame(value=d[rownwames(d)=="DEFA5",])
pd <- data.frame(value=d[rownames(d)=="DEFA5",])
ggplot(data=pd, aes(value)) +  geom_histogram(fill="darkgreen")+theme_bw()
heaD(pd)
head(pd)
pd <- data.frame(value=as.numeric(d[rownames(d)=="DEFA5",])
)
head(pd)
ggplot(data=pd, aes(value)) +  geom_histogram(fill="darkgreen")+theme_bw()
min(d)
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4')
#load('')
load('dds.Rdata')
rownames(pd) <- colnames(d)
head(pd)
pd <- pd[grepl('.',rownames(pd), fixed=TRUE),]
head(pd)
pd <- pd[grepl('.',rownames(pd), fixed=TRUE),, drop=F]
pd <- pd[grepl('.',rownames(pd), fixed=TRUE), drop=F]
head(pd)
pd <- data.frame(value=as.numeric(d[rownames(d)=="DEFA5",]))
rownames(pd) <- colnames(d)
table(grepl('.',rownames(pd), fixed=TRUE))
pd <- pd[!grepl('.',rownames(pd), fixed=TRUE),, drop=F]
head(pd)
head(fpkm_d)
head(pd)
d <- data.frame(value=fpkm_d[rownames(fpkm_d)=="H_DEFA5",])
rownames(d) <- colnames(fpkm_d)
d <- d[grepl("LMX", rownames(d)), ,drop=F]
head(d)
rownames(d) <- substr(rownames(d), 0,9)
d$case  <- substr(rownames(d), 0,9)
table(d$case %in% rownames(metadata[metadata$type=="BASALE",])))
table(d$case %in% rownames(metadata[metadata$type=="BASALE",]))
table(d$case %in% rownames(metadata[grepl("BASALE",metadata$type),]))
head(metadta)
head(metadata)
table(d$case %in% rownames(metadata[grepl("BASALE",metadata$type, perl=TRUE),]))
table(d$case %in% rownames(metadata[grepl("BASALE",metadata$type, perl=TRUE),]))
table(grepl("BASALE",metadata$type, perl=TRUE))
head(rownames(metadata[grepl("BASALE",metadata$type, perl=TRUE))
head(rownames(metadata[grepl("BASALE",metadata$type, perl=TRUE),])
)
head(d)
table(rownames(d) %in% rownames(metadata[grepl("BASALE",metadata$type),]))
d <- d[rownames(d) %in% rownames(metadata[grepl("BASALE",metadata$type),])),]
d <- d[rownames(d) %in% rownames(metadata[grepl("BASALE",metadata$type),]),]
head(d)
head(pd)
m <- merge(d, pd, by.x="case", by.y="row.names")
head(m)
cor.test(value.x, value.y)
cor.test(m$value.x, m$value.y)
plot(m$value.x, m$value.y)+xlab('RNAseq')+ylab('uarray')
plot(m$value.x, m$value.y,xlab='RNAseq',ylab='uarray')
m$value.x <- log(m$value.x+0.0001)
plot(m$value.x, m$value.y,xlab='RNAseq',ylab='uarray')
cor.test(m$value.x, m$value.y,xlab='RNAseq',ylab='uarray')
m[m$value.x < -5,]
d
d$l <- log(d$value+0.0001)
d[d$l < -5]
d[d$l < -5,]
###
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4')
#load('')
load('dds.Rdata')
pdxbasali <- grepl("BASALE",metadata$type, perl=TRUE) & grepl('LMX', rownames(metadata), perl=TRUE)
table(pdxbasali)
f <- fpkm[,pdxbasali]
f <- fpkm_d[,pdxbasali]
#substr(0,9)
grepl('.',colnames(d), fixed=TRUE)
rep <- grepl('.',colnames(d), fixed=TRUE)
d <- d[,!res]
d <- as.data.frame(d)
d <- d[,!res]
d <- read.table(gzfile('~/matrice_microarray.tsv.gz'), sep="\t", header=T)
head(d)
rep <- grepl('.',colnames(d), fixed=TRUE)
head(rep)
d <- d[,!rep]
head(d)
keep <- intersect(colnames(d), substr(colnames(f,0,9)))
keep <- intersect(colnames(d), substr(colnames(f),0,9))
length(keep)
colnames(f) <- substr(colnames(f),0,9)
f <- f[, colnames(f) %in% keep]
d <- d[, colnames(d) %in% keep]
dim(d)
dim(f)
)
rownames(f) <- substr(rownames(f),2,nchar(rownames(f)))
keepr <- intersect(rownames(f), rownames(d))
length(keepr)
head(rownames(f))
rownames(f) <- substr(rownames(f),1,nchar(rownames(f)))
head(rownames(f))
rownames(f) <- substr(rownames(f),2,nchar(rownames(f)))
head(rownames(f))
keepr <- intersect(rownames(f), rownames(d))
length(keepr)
d <- d[rownames(d) %in% keepr,]
f <- f[rownames(f) %in% keepr,]
head(f[rownames(d),])
head(rownames(d))
f <- f[rownames(d),])
f <- f[rownames(d),]
head(f[,colnames(d)])
head(d)
f <- f[,colnames(d)]
all(rownames(f)==rownames(d))
all(colnames(f)==colnames(d))
dim(d)
f <- log(f+0.0001)
#correl <- sapply(colnames(f), function(f) {cor()})
correl <- cor(f, d)
head(correl)
plot(correl)
#correl <- sapply(rownames(f), function(f) {cor()})
correl <- sapply(rownames(f), function(f) {cor()})
correl <- cor(t(f), t(d))
head(correl)
correl[names(correl)=="DEFA5"]
head(names(correl))
correl[rownames(correl)=="DEFA5"]
correl[rownames(correl)=="DEFA5",]
dim(correl)
correl[rownames(correl)=="DEFA5", colnames(correl)=="DEFA5"]
plot(d[rownames(d)=="DEFA5"], f[rownames(f)=="DEFA5"])
plot(d[rownames(d)=="DEFA5",], f[rownames(f)=="DEFA5",])
hea(d)
head(d)
x1 <- d[rownames(d)=="DEFA5",]
y1 <- f[rownames(f)=="DEFA5",]
head(x1)
plot(as.numeric(x1), as.numeric(y1))
plot(as.numeric(y1), as.numeric(x1))
cor.test(as.numeric(y1), as.numeric(x1))
dim(correl)
diagc <- diag(correl)
length(diagc)
plot(density(diagc))
head(diagc))
head(diagc)
plot(density(as.numeric(diagc)))
summary(diagc)
diagc <-diagc[!is.na(diagc)]
summary(diagc)
plot(density(as.numeric(diagc)))
head(diagc[order(diagc)])
head(diagc[order(-diagc)])
plot(as.numeric(f[rownames(f)=="ARHGAP4",]), as.numeric(d[rownames(d)=='ARHGAP4',]))
plot(as.numeric(f[rownames(f)=="IL10",]), as.numeric(d[rownames(d)=='IL10',]))
summary(diagc)
plot(as.numeric(f[rownames(f)=="MIR2116",]), as.numeric(d[rownames(d)=='MIR2116',]))
### paired de
load('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4_cetuxi_treat_PDX/dds.Rdata')
library(DESeq2)
head(colnames(data))
head(rownames(metadata))
head(colnames(newdata))
head(colnames(new_data))
metadata$tt <- paste0(metadata$time, metadata$treat)
fdesign <- as.formula('batch+sample+tt')
fdesign <- as.formula('~batch+sample+tt')
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
head(metadata)
metadata$batch <- as.factor(metadata$batch)
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
fdesign
fdesign <- as.formula("~batch+sample+time+treat")
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
fdesign <- as.formula("~sample+time+treat")
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
metadata
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatf <- paste0(outprefix, "_vsd_dists.pdf")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, file=pheatf)
pheatf
heatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off*()
dev.off()
dev.off()
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, annotation_row = metadata[,c("time","treat")])
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, annotation_row = metadata[,c("time","treat","sample")])
plotPCA(vsd, intgroup='treat')
plotPCA(vsd, intgroup='treat')
ints <- attr(terms(fdesign),"term.labels")
ints
plotPCA(vsd, intgroup='sample)
'
)
plotPCA(vsd, intgroup='sample')
head(vsd)
library(DESeq2)
head(vsd)
plotPCA(vsd, intgroup='sample')
plotPCA
plotPCA(dds, intgroup='sample')
exprs(dds)
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
ppca(vsd, 'treat',1,2)
ppca(vsd, 'time',1,2)
ppca(vsd, 'sample',1,2)
dds2 <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = as.formula('~batch'))
filterGenes <- rowSums(counts(dds2) > minc) < minsamples
dds2 <- dds2[!filterGenes]
dds2 <- DESeq(dds2, parallel=TRUE, betaPrior=TRUE)
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)
vsd2 <- vst(dds2, blind=FALSE)
ppca(vsd, 'time',1,2)
ppca(vsd2, 'time',1,2)
ppca(vsd2, 'treat',1,2)
ppca(vsd, 'treat',1,2)
ppca(vsd, 'batch',1,2)
ppca(vsd2, 'batch',1,2)
ppca(vsd2, 'sample',1,2)
ppca(vsd, 'sample',1,2)
alpha <- 0.05
lfc <- 0.5849625
con <- c('treat', 'cetuxi','NT')
parallel <- TRUE
res <- results(dds, alpha=alpha, contrast=con, parallel=parallel)
res
plot_volcano <- function(resnona, alpha, lfc, outfile, title) {
#title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
#resnona <- res[!is.na(res$padj),]
resnona$sign <- ifelse(abs(resnona$log2FoldChange) > lfc & resnona$padj < alpha, "both", ifelse(abs(resnona$log2FoldChange) > lfc, "LFC",
ifelse(resnona$padj < alpha, "padj", "NS")))
resnona$sign <- factor(resnona$sig, levels=c("LFC", "padj", "both", "NS"))
resnona[resnona$padj ==0,"padj"] <- .Machine$double.xmin
p <- ggplot(resnona, aes(log2FoldChange, -log10(padj))) +
geom_point(aes(col = sign),size=0.5) + theme_bw() +
scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE) + # red orange green black -> orange blue green gray
ggtitle(title)
nsign <- nrow(resnona[resnona$sig=="both",])
if (nsign > 20) {
p + geom_text_repel(data=resnona[1:10,], aes(label=rownames(resnona)[1:10]))
} else {
p + geom_text_repel(data=resnona[resnona$sig=="both",], aes(label=rownames(resnona[resnona$sig=="both",])))
}}
resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])
title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
plot_volcano(resnona_df, alpha, lfc, volcano, title)
write.table(resnona_df, file='de_paired.tsv', quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
res2 <- read.table('treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', sep="\t", header=T)
dim(res)
dim(resnona_df)
dim(res2)
m <- merge(res2, resnona_df, by="row.names")
head(m)
plot(m$baseMean.x,m$baseMean.y)
plot(m$log2FoldChange.y,m$log2FoldChange.x)
plot(m$padj.y,m$padj.x)
m1 <- m[m$padj.x < 0.01,]
m2 <- m[m$padj.y < 0.01,]
library(VennDiagram)
venn.diagram(
x = list(rownames(m2), rownames(m1)),
category.names = c("unpaired" , "paired"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
x = list(rownames(m2), rownames(m1)),
category.names = c("unpaired" , "paired"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
plot(log10(m$padj.y),log10(m$padj.x))
plot(log10(m$padj.y),log10(m$padj.x), xlabs="unpaired",ylabs="paired")
plot(log10(m$padj.y),log10(m$padj.x), xlabs="unpaired",ylabs="paired")
plot(log10(m$padj.y),log10(m$padj.x), xlab="paired",ylab="unpaired")
warnings()
venn.diagram(
x = list(rownames(m2), rownames(m1)),
category.names = c("paired" , "unpaired"),
output=TRUE,
# Output features
imagetype="png" ,
filename="venn_lfc.png", fill = myCol[c(1,2)]
)
m[rownames(m)=='H_DEFA5',]
m[rownames(m)=='DEFA5',]
head(m)
m[m$Row.names=='H_DEFA5',]
m[m$Row.names=='H_DEFA6',]
wanted <- c("H_DEFA5","H_DEFA6",'H_ATOH1','H_DLL1','H_GFI1','H_SOX9','H_SPINK4','H_XBP1','H_DLL4','H_LYZ')
wanted2 <- c("H_TERT","H_HOPX",'H_BMI1','H_LRIG1')
library(ggplot2)
mm <- m[m$Row.names %in% wanted,]
mm
ggplot(data=mm, aes(x=reorder(Row.names, -log2FoldChange.y), fill=log10(padj.y),y=log2FoldChange.y)) +
geom_bar(stat="identity")+theme_bw()+scale_fill_distiller(palette="Spectral")+ggtitle('Paneth')+theme(text = element_text(size=20))

mm2 <- m[m$Row.names %in% wanted2,]
ggplot(data=mm2, aes(x=as.factor(Row.names), fill=log10(padj.y),y=log2FoldChange.y)) +
geom_bar(stat="identity")+theme_bw()+scale_fill_distiller(palette="Spectral")+ggtitle('LRC')+theme(text = element_text(size=20))


