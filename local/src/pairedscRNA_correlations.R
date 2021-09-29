d1 <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/scRNA_paired/vsd.tsv.gz'), sep="\t")
d2 <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz'), sep="\t")
meda <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data'), sep="\t", header=TRUE)
meda <- meda[grepl("LMO_cetuxi", meda$type), ]
d3 <- d2[,colnames(d2) %in% meda$id]
rownames(d3) <- substr(rownames(d2), 3, nchar(rownames(d2)))
md <- merge(d1, d3, by="row.names")
rownames(md) <- md$Row.names
md$Row.names <- NULL

means <- apply(md, 1, mean)
med <- median(means)
# 
de <- md[means > med,]
 
sds <- apply(de, 1, sd)
# # there are some very high sd?
# #    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# #0.9      5.6     28.9    354.0    149.4 860110.9 
# noisy_genes_thr <- quantile(sds, 0.9)
# 
# de <- de[sds < noisy_genes_thr,]
# sds <- sds[sds < noisy_genes_thr]
# 
# # now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]
# 
n <- length(sds)
keep <- head(sds, round(0.10*n))
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

meda2 <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/scRNA_paired/samples_data'), sep="\t", header=TRUE)

meda$cetuxi <- NULL
meda$irino <- NULL
meda$batch <- NULL
meda$treat <- ifelse(grepl("NT", meda$type), 'NT_72h','Cetux_72h')
meda$sort <- 'pre-sorting'
meda$type <- NULL
meda$smodel <- substr(meda$id, 0, 7)

ameda <- rbind(meda, meda2)
rownames(ameda) <- ameda$id
ameda$id <- NULL


library(pheatmap)

pheatmap(as.matrix(desd), annotation_col = ameda)

keep <- c('CRC0542','CRC0069', 'CRC1307')

ameda2 <- ameda[ameda$smodel %in% keep,]
desd2 <- desd[, colnames(desd) %in% rownames(ameda2)]
pheatmap(as.matrix(t(desd2)), annotation_row = ameda2, show_colnames = F)
