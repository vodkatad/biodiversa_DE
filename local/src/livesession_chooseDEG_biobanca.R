load('/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584/dds.Rdata')
#/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/trialclustering/LMO.vs.LMH_allgenes_category.tsv 

data <- plotCounts(dds, gene="H_CHEK2", intgroup="type", returnData=TRUE)
data$sample <- substr(rownames(data), 0, 7)
ggplot(data=data, aes(x=type, y=count, color=sample))+geom_point()+scale_y_log10()+theme_bw()+geom_line(aes(group=sample))+ theme(legend.position = "none")



data <- plotCounts(dds, gene="H_PDX1", intgroup="type")
data <- plotCounts(dds, gene="H_PDX1", intgroup="type", returnData=TRUE)
data$sample <- substr(rownames(data), 0, 7)
ggplot(data=data, aes(x=type, y=count, color=sample))+geom_point()+scale_y_log10()+theme_bw()+geom_line(aes(group=sample))+ theme(legend.position = "none")


gene <- 'H_DEFA6'
plotCounts(dds, gene=gene, intgroup="type")
data <- plotCounts(dds, gene=gene, intgroup="type", returnData=TRUE)
data$sample <- substr(rownames(data), 0, 7)
ggplot(data=data, aes(x=type, y=count, color=sample))+geom_point()+scale_y_log10()+theme_bw()+geom_line(aes(group=sample))+ theme(legend.position = "none")

gene <- 'H_SC5D'
plotCounts(dds, gene=gene, intgroup="type")
data <- plotCounts(dds, gene=gene, intgroup="type", returnData=TRUE)
data$sample <- substr(rownames(data), 0, 7)
ggplot(data=data, aes(x=type, y=count, color=sample))+geom_point()+scale_y_log10()+theme_bw()+geom_line(aes(group=sample))+ theme(legend.position = "none")


gene <- 'H_REN'
plotCounts(dds, gene=gene, intgroup="type")
data <- plotCounts(dds, gene=gene, intgroup="type", returnData=TRUE)
data$sample <- substr(rownames(data), 0, 7)
ggplot(data=data, aes(x=type, y=count, color=sample))+geom_point()+scale_y_log10()+theme_bw()+geom_line(aes(group=sample))+ theme(legend.position = "none")


### replicates correlations
datax <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_replicates_correlation.tsv', sep="\t", header=T)
datao <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_replicates_correlation.tsv', sep="\t", header=T)
datax$class <- 'LMX'
datao$class <- 'LMO'
data <- rbind(datax, datao)
library(ggplot2)
ggplot(data=data, aes(x=correlation, fill=class))+geom_histogram(alpha=0.5)+theme_bw()+theme(text=element_text(size=20))
#ggplot(data=data, aes(x=correlation, fill=class))+geom_histogram(alpha=0.5)+theme_bw()+theme(text=element_text(tsize=20))
table(data$class)

### stromal investigations
d <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/PDO_buoni_Cutoff0.05_LFC0.584/type_cutoff0.05-LMX_BASALE.vs.LMH.deseq2.tsv', sep="\t", header=TRUE)
s <- read.table('/scratch/trcanmed/DE_RNASeq/local/share/data/Stromal_genes.tsv', sep="\t")
d$gene <- gsub('H_', '', rownames(d))
d$isstromal <- d$gene %in% s$V1
table(d$isstromal)

ggplot(data=d, aes(y=log2FoldChange, x=isstromal, fill=isstromal))+geom_boxplot()+theme_bw()+theme(text=element_text(size=20))

wilcox.test(d[!d$isstromal, 'log2FoldChange'],d[d$isstromal, 'log2FoldChange'])


  d <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/PDO_buoni_Cutoff0.05_LFC0.584/type_cutoff0.05-LMO_BASALE.vs.LMH.deseq2.tsv', sep="\t", header=TRUE)
  d$gene <- gsub('H_', '', rownames(d))
  d$isstromal <- d$gene %in% s$V1
  table(d$isstromal)
  
  
load('/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584/dds.Rdata')
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
counts <- counts(dds, normalized=TRUE)

#(bit_rnaseq_2.8) egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/PDO_buoni_Cutoff0.05_LFC0.584$ filter_1col 1  <(cut -f 1 /scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/gene_len| sed 's/H_//1' ) < /tmp/stromali_nouniverso  > /tmp/daspiegare
dd <- read.table('/tmp/daspiegare')

#d <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/PDO_buoni_Cutoff0.05_LFC0.584/type_cutoff0.05-LMX_BASALE.vs.LMH.deseq2.tsv', sep="\t", header=TRUE)
#s <- read.table('/scratch/trcanmed/DE_RNASeq/local/share/data/Stromal_genes.tsv', sep="\t")
#d$gene <- gsub('H_', '', rownames(d))
#d$isstromal <- d$gene %in% s$V1
#d$isstromal_have <- d$gene %in% dd$V1

rownames(counts) <- gsub('H_', '', rownames(counts))
cstro <- counts[rownames(counts) %in% dd$V1,]

library(pheatmap)
annot <- metadata[, 'type', drop=FALSE]
annot <- annot[order(annot$type), , drop=FALSE]
cstro <- cstro[,rownames(annot)]
pheatmap(log2(cstro+1), annotation_col=annot, show_rownames = F, show_colnames=F, cluster_cols=F)


