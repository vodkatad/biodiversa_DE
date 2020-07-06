d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4/all_prediction_result.xls', sep="\t", header=T)
w <- read.table('~/CHEMIO_WATERFALL_PLOT_Eugy_2020maggio_w3.txt', sep="\t", header=T)
head(w)
head(d)
dd <- d[grepl('LMX',d$sample.names,fixed=T),]
head(dd)
dd$case <- substr(dd$sample.names, 0, 7)
head(dd)
head(w)
ddd <- dd[dd$BH.FDR < 0.05,]
dim(ddd)
dim(dd)
m <- merge(ddd, w, by="case")
head(m)
ggplot(m, aes(x=predict.labal2, y=perc)) +
geom_boxplot()
ggplot(m, aes(x=predict.label2, y=perc)) +
geom_boxplot()


> ss <- substr(basali$id, 0, 7)
> length(ss)
[1] 266
> head(ss)
[1] "CRC0018" "CRC0024" "CRC0024" "CRC0028" "CRC0029" "CRC0050"
> length(unique(ss))
[1] 236


d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4/all_prediction_result.xls', sep="\t", header=T)
w <- read.table('~/CHEMIO_WATERFALL_PLOT_Eugy_2020maggio_w3.txt', sep="\t", header=T)
me <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4/samples_data', sep="\t", header=T)
basali <- me[grepl('LMX_BASALE', me$type),]
dim(basali)
head(basali)
basali
ss <- substr(basali$id, 0, 7)
length(ss)
head(ss)
length(unique(ss))
ddd <- dd[dd$BH.FDR < 0.2,]
ddd <- d[d$BH.FDR < 0.2,]
head(ddd)
head(w)
ddd$case <- substr(ddd$sample.names, 0, 7)
length(ddd$cse)
length(ddd$case)
basali_cris <- merge(ddd, basali, by.x='sample.names',by.y='id')
head(basali_cris)
dir(basali_cris)
dim(basali_cris)
length(basali_cris$type)
length(unique(basali_cris$type))
length(unique(basali_cris$case))
length(basali_cris$case)
table(basali_cris$case)
dups <- table(basali_cris$case)
dups_cases <- names(dups[dups==2])
dups_cases
basali_cris[basali_cris$case %in% dups_cases, ]
dups_cris <- sapply(dups_cases, function(x) { y <- basali_cris[basali_cris$case %in% x }; length(unique(y$predict.label2))]
dups_cris <- sapply(dups_cases, function(x) { y <- basali_cris[basali_cris$case %in% x,]; length(unique(y$predict.label2)) }
)
dups_cris
remove <- names(dups_cris[dups_cris==2,])
remove <- names(dups_cris[dups_cris==2])
remove
basali_cris_nodup <- basali_cris[!basali_cris$case %in% remove,]
dim(basali_cris_nodup)
dim(basali_cris)
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris[basali_cris$case %in% x,]; y$sample_names[1] }
)
same_cris
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris[basali_cris$case %in% x,]; y$sample_names }
)
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris_nodup[basali_cris_nodup$case %in% x,]; y$sample_names } )
head(same_cris)
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris_nodup[basali_cris_nodup$case %in% x,]; y$sample.names } )
head(same_cris)
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris_nodup[basali_cris_nodup$case %in% x,]; y$sample.names[1] } )
head(same_cris)
ssame_cris <- same_cris[same_cris != NA]
head(ssame_cris)
same_cris <- sapply(dups_cases, function(x) { y <- basali_cris_nodup[basali_cris_nodup$case %in% x,]; as.character(y$sample.names[1]) } )
head(same_cris)
ssame_cris <- same_cris[!is.na(same_cris)]
ssame_cris
basali_cris_nodup_onerep <- basali_cris_norup[basali_cris_nodup$sample.names %in% ssame_cris,]
basali_cris_nodup_onerep <- basali_cris_nodup[basali_cris_nodup$sample.names %in% ssame_cris,]
dim(basali_cris_nodup_onerep)
basali_cris_nodup_onerep <- basali_cris_nodup[!basali_cris_nodup$sample.names %in% ssame_cris,]
dim(basali_cris_nodup_onerep)
ggplot(m, aes(x=predict.labal2, y=cetuxi)) +
geom_boxplot()
ggplot(basali_cris_nodup_onerep, aes(x=predict.labal2, y=cetuxi)) +
geom_boxplot()
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=cetuxi)) +
geom_boxplot()
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=cetuxi, fill=predict.label2)) +    geom_boxplot()+theme_bw()
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=cetuxi, fill=predict.label2)) +    geom_boxplot()+theme_bw()+ggtitle('cetuxi w3')
head(basali_cris_nodup_onerep)
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=cetuxi*100, fill=predict.label2)) +    geom_boxplot()+theme_bw()+ggtitle('cetuxi w3')
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=irino*100, fill=predict.label2)) +    geom_boxplot()+theme_bw()+ggtitle('cetuxi w3')
m <- merge(basali_cris_nodedup_onerep, w, by="case")
m <- merge(basali_cris_nodup_onerep, w, by="case")
head(m)
dim(m)
dim(basali_cris_nodup_onerep[!is.na(basali_cris_nodup_onerep$irino),])
dim(basali_cris_nodup_onerep[!is.na(basali_cris_nodup_onerep$cetuxi),])
ggplot(m, aes(x=predict.label2, y=perc, fill=predict.label2)) +    geom_boxplot()+theme_bw()+ggtitle('irino w3')
ggplot(m, aes(x=predict.label2, y=perc, fill=predict.label2)) +    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2)+theme_bw()+ggtitle('irino w3')+
)
ggplot(m, aes(x=predict.label2, y=perc, fill=predict.label2)) +geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2)+theme_bw()+ggtitle('irino w3')
ggplot(basali_cris_nodup_onerep, aes(x=predict.label2, y=cetuxi*100, fill=predict.label2)) +geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2)+theme_bw()+ggtitle('cetuxi w3')

#### lympho scores
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4')
#load('')
load('dds.Rdata')
CA(vsd, type, returnData=T)
pcs <- DESeq2::plotPCA(vsd, 'type', returnData=T)

pc2_42 <- pcs[pcs$PC2>42,]

stromal <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4/stromal_contamination.txt', header=T, sep="\t")
stromal <- stromal[, c(1,2,3,4)]
m <- merge(stromal, pcs, by.x='X', by.y="row.names")
plot(m$leucos, m$PC2)
plot(m$cafs, m$PC2)
plot(m$endos, m$PC2)
quantile(stromal$leucos, probs=0.9)
dim(stromal[stromal$leucos>3.58,])

lymph_large <- stromal[stromal$leucos>2.23,]
lymph_strict <- stromal[stromal$leucos>5,]

lymph_strict$case <- substr(as.character(lymph_strict$X), 0,10)
pc2_42$case <- substr(as.character(rownames(pc2_42)), 0,10)
d <- data.frame(ids=c(lymph_strict$case,pc2_42$case))
write.table(d, file="lympho_pc_stromal_strict", sep="\t", col.names = F, quote=F, row.names=F)
write.table(pc2_42, file="pc2_42.tsv", row.names=T, col.names=T, quote=F, sep="\t")


### normalized counts
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up4')
#load('')
load('dds.Rdata')
h_counts <- counts(dds, normalized=TRUE)
data <- read.table(gzfile(counts), header=T, sep="\t", row=1)
prefix <- "M"
metadata <- read.table(metadataf, sep="\t", header=T, row=1)
fdesign <- as.formula(design)
print(terms(fdesign)[[2]])
new_data <- data[,match(rownames(metadata), colnames(data))]
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)
m_counts <- counts(dds, normalized=TRUE)
head(m_counts)
dim(h_counts)
dim(m_counts)
save.image('counts_h_m_corrected.Rdata')
plotPCA(vsd, intgroup='type')
DESeq2::plotPCA(vsd, intgroup='type')
design
#common <- intersect(rownames(h_counts), rownames(m_counts))
counts <- rbind(h_counts, m_counts)
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
mousemean <- colMeans(m_mouse)
mousemean <- colMeans(m_counts)
df <- data.frame(rownames=names(mousemean), count=mousemean)
head(df)
df$type <- substr(rownames(df), 8,10)
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+scale_y_log10()
library(ggsignif)
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMX'), c("LMO","LMX"),c('PRH','PRX')))+scale_y_log10()
df$hx <- ifelse(df$type %in% c("LMX",'PMX','PRX'), 'X', 'H')
ggplot(df, aes(x=hx, y=count, fill=hx)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("X","H"))
)+scale_y_log10()
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMX'), c("LMO","LMX"),c('PRH','PRX')))
ggplot(df, aes(x=type, y=count, fill=type)) + geom_boxplot()+theme_bw()+geom_signif( map_signif_level=TRUE, comparisons=list(c("LMH",'LMX'), c("LMO","LMX"),c('PRH','PRX')))+ylim(0,1000)


