library(DESeq2)
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S')
load('dds.Rdata')
nt <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/pb/VandE_CRC0327_NT_2_bulklog2.csv.tsv', sep='\t', header=F, row.names=1)
ct <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/pb/VandE_CRC0327_cetux_2_bulklog2.csv.tsv', sep='\t', header=F, row.names=1)

nt2 <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/pb/VandACRC0322_NT_1_bulklog2.csv.tsv', sep='\t', header=F, row.names=1)
ct2 <- read.table('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/pb/VandECRC0322cetux1_bulklog2.csv.tsv', sep='\t', header=F, row.names=1)

scdata2 <- data.frame(CRC0327_ctx2=round(rowSums(exp(ct)+1)), CRC0327_nt2=round(rowSums(exp(nt)+1)), CRC0322_ctx1=round(rowSums(exp(ct2)+1)), CRC0322_nt=round(rowSums(exp(nt2)+1)))
rownames(scdata2) <- paste0("H_", rownames(scdata2))
m <- merge(data, scdata2, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL
metadata <- read.table(metadataf, sep="\t", header=T, row=1)
metadata <- rbind(metadata, data.frame(row.names=c('CRC0327_ctx2','CRC0327_nt2','CRC0322_ctx1','CRC0322_nt1'), sample=c('CRC0327','CRC0327','CRC0322','CRC0322'),batch=c(6,6,6,6), treat=c('cetuxi','NT','cetuxi','nt')))
#metadata <- rbind(metadata, data.frame(row.names=c('CRC0327_ctx2','CRC0327_nt2'), sample=c('CRC0327','CRC0327'),batch=c(6,6), treat=c('cetuxi','NT')))
metadata[,'batch'] <- as.factor(metadata[,'batch'])

new_data <- m[,match(rownames(metadata), colnames(m))]
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup='treat')
plotPCA(vsd, intgroup='sample')
plotPCA(vsd, intgroup='batch')
metadata <- metadata[metadata$sample=="CRC0327",]
new_data <- m[,match(rownames(metadata), colnames(m))]

metadata$sample <- NULL
fdesign <- as.formula("~treat+batch")
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)



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


lens <- read.table(len, sep="\t", header=TRUE)
order <- rownames(mcols(dds))
lens <- lens[match(order, lens$Geneid), ]
mcols(dds)$basepairs  <- lens$length
fpkm_d <- fpkm(dds)

res <- results(dds, alpha=0.05, contrast=c('treat','cetuxi','NT'))
resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])

paneth <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/paneth_s_r_signature'), sep="\t", header=F)
colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)

f <- fpkm_d[rownames(fpkm_d) %in% paneth$hgs,]
pheatmap(log(f+0.0001), annotation_col = metadata)

fcs <- data.frame(fc_bulk1=fpkm_d[,2]/fpkm_d[,1], fc_bulk2=fpkm_d[,4]/fpkm_d[,3], fc_sc=fpkm_d[,5]/fpkm_d[,6])
fcs <- data.frame(fc_bulk1=fpkm_d[,2]/(fpkm_d[,1]+0.0001), fc_bulk2=fpkm_d[,4]/(fpkm_d[,3]+0.0001), fc_sc=fpkm_d[,5]/fpkm_d[,6])
fcss <- fcs[rownames(fcs) %in% paneth$hgs,]
pheatmap(log2(fcss+1), show_rownames=TRUE, show_colnames=TRUE)
sign <- resnona[resnona$padj< 0.05 & resnona$log2FoldChange>1.5,]
dim(sign)
fcss <- fcs[rownames(fcs) %in% rownames(sign),]
pheatmap(log2(fcss+0.0001), show_rownames=TRUE, show_colnames=TRUE)
sign <- resnona[resnona$padj< 0.05 & resnona$log2FoldChange< -1.5,]
dim(sign)
sign <- resnona[resnona$padj< 0.05 & resnona$log2FoldChange< -2,]
dim(sign)
fcss <- fcs[rownames(fcs) %in% rownames(sign),]
pheatmap(log2(fcss+0.0001), show_rownames=TRUE, show_colnames=TRUE)
summary(fcss)
summary(fcs)

v <- assay(vsd)
tc <- cor(v)
corrplot(tc)

###


res2 <- results(dds, alpha=0.05, contrast=c('treat','cetuxi','NT'))
resnona2 <- res2[!is.na(res2$pvalue) & !is.na(res2$padj),]
resnona_df2 <- as.data.frame(resnona2[order(resnona2$padj),])