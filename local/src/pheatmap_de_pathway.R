setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Dambrosio_LMX_DE')
library('pheatmap')
#fpkm <- read.table(gzfile('fpkm.tsv.gz'), header=T, sep="\t")
tmm <- read.table(gzfile('tmm.tsv.gz'), header=T, sep="\t")
deg <- read.table('mut_cutoff0.05-MUT.vs.WT.deseq2_sign_genedesc.tsv', header=T, sep="\t", quote='', comment.char='')
tnf <- read.table('../../local/share/data/HALLMARK_TNFA_SIGNALING_VIA_NFKB', header=F)
#tnf <- read.table('../../local/share/data/REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE', header=F) # seems true even if less expressed
#tnf <- read.table('../../local/share/data/REACTOME_RECRUITMENT_OF_NUMA_TO_MITOTIC_CENTROSOMES', header=F) # ok does not seem true
tnf_other <- read.table('../../local/share/data/REACTOME_TNF_SIGNALING', header=F)
length(intersect(deg$gene, tnf$V1))
length(intersect(deg$gene, tnf_other$V1))
intersect(deg$gene, tnf_other$V1)
intersect(deg$gene, tnf$V1)
length(intersect(tnf_other$V1, tnf$V1))
length(tnf_other$V1); length(tnf$V1)
intersect(tnf_other$V1, tnf$V1)
means <- rowMeans(tmm)
pdata <- data.frame(rownames=rownames(tmm), means=means)
pdata$meandeciles <- cut( pdata$mean, quantile(pdata$mean, prob = seq(0, 1, length = 11), type = 5), include.lowest=TRUE )
levels(pdata$meandeciles) <- seq(0,1, length=11)

summary(pdata[rownames(pdata) %in% paste0("H_",tnf$V1),'meandeciles'])

rownames(tmm) <- substr(rownames(tmm),3, nchar(rownames(tmm)))
#rownames(fpkm) <- substr(rownames(fpkm),3, nchar(rownames(fpkm)))
anno <- read.table('samples_data', sep="\t", header=T)
annot_col <- data.frame(row.names=anno$id, KRAS=anno$mut)

pdata <- tmm[rownames(tmm) %in% tnf$V1,]

annot_rows <- data.frame(row.names=rownames(pdata))
annot_rows$DEG <- ifelse(rownames(annot_rows) %in% deg$gene, 'yes','no')

#pheatmap(log10(pdata+0.0001), annotation_row = annot_rows, annotation_col=annot_col)

#pzscore <- zscore[rownames(zscore) %in% tnf$V1,]
means <- apply(tmm, 1, mean)
tmm <- tmm[!rownames(tmm) %in% names(means[means==0]),]
ltmm <- log(tmm+0.0001)
z <- t(scale(t(ltmm)))
pz <- z[rownames(z) %in% tnf$V1,]
pheatmap(pz, annotation_row = annot_rows, annotation_col=annot_col)

pdata <- tmm[rownames(tmm) %in% tnf$V1,]

#pheatmap(log10(pdata+0.0001), annotation_row = annot_rows, annotation_col=annot_col)

all <- read.table('mut_cutoff0.05-MUT.vs.WT.deseq2.tsv', sep="\t", header=T)
rownames(all) <- substr(rownames(all),3, nchar(rownames(all)))
lll <- all[rownames(all) %in% tnf$V1,]
lll$padj <- p.adjust(lll$pvalue, method="BH")
lll <- lll[lll$padj < 0.05,]
selected <- rownames(lll)

pdata <- pdata[rownames(pdata) %in% selected, ]
annot_col <- annot_col[order(annot_col$KRAS),, drop=F]
pdata <- pdata[, match(rownames(annot_col), colnames(pdata))]

pheatmap(log10(pdata+0.0001), annotation_col=annot_col, cluster_cols = F, show_colnames=F)


pdata <- pz[rownames(pz) %in% selected, ]
pdata <- pdata[, match(rownames(annot_col), colnames(pdata))]
pheatmap(pdata, annotation_col=annot_col, cluster_cols = F, show_colnames=F)

### deg boxplots
pdata <- tmm[rownames(tmm) %in% deg$gene,]
mm <- median(deg$baseMean)
degg <- deg[deg$log2FoldChange > 1.5,]
pdata <- tmm[rownames(tmm) %in% degg$gene,]
#pheatmap(log10(pdata+0.0001), annotation_col=annot_col, cluster_cols = T, show_colnames=F, show_rownames=F)

#
pdata <- tmm[rownames(tmm) %in% tnf$V1,]
pdata <- pdata[rownames(pdata) %in% selected, ]
td <- t(pdata)
long <- melt(td)
mlong$x <- paste0( mlong$X2, "_", mlong$KRAS)
all$gene <- rownames(all)
d <- all[all$gene %in% long$X2,]
#ggplot(data=mlong, aes(x=x,y=value,color=X2))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1))

order <- d[order(-d$log2FoldChange),]
f <- paste0(order$gene, "_", 'MUT')
f2 <- paste0(order$gene, "_", 'WT')
ff <- c()
for(i in seq(1, length(f2))) {
  ff <- c(ff, f[i])
  ff <- c(ff, f2[i])
}
mlong$x <- factor(mlong$x, levels=ff)
ggplot(data=mlong, aes(x=x,y=value,color=X2))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1))

ggplot(data=mlong, aes(x=x,y=value,color=X2))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1))+scale_y_log10()