meta <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/general/samples_data', sep="\t", header=T, stringsAsFactors = F)


d <- read.table(gzfile('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/general/tmm.tsv.gz'), sep="\t", header=T)

#paneth <- data.frame(gs=c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4"))
paneth <- data.frame(gs=c("ATOH1","GFI1","DEFA5","DEFA6","DLL1"))

colnames(paneth) <- 'gs'

fpkm <- d[rownames(d) %in% paneth$gs,]
pseudoc <- 1
lfpkm <- log2(fpkm+pseudoc)

library(pheatmap)
pheatmap(lfpkm, annotation_col = meta)


ave <- colMeans(lfpkm)
paneth_treat <- data.frame(ave =ave, id = names(ave), stringsAsFactors = FALSE)
paneth_treat_m <- paneth_treat[grepl('WT', paneth_treat$id),]
paneth_treat_m$model <- substr(paneth_treat_m$id, 0, 7)
fc <- function(model, data, col1, col2) {
  d <- data[data$model == model,]
  fc <- d[grepl(col1, d$id), 'ave'] - d[grepl(col2, d$id), 'ave']
  return(fc)
}

panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'CTX', 'NO_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))




panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'NO_EGF', 'WT_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))




panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'CTX', 'WT_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))

############ PIS CAS9 TODO


# DLL1 only
meta <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/general/samples_data', sep="\t", header=T, stringsAsFactors = F)


d <- read.table(gzfile('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/general/tmm.tsv.gz'), sep="\t", header=T)

#paneth <- data.frame(gs=c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4"))
paneth <- data.frame(gs=c("DLL1"))

colnames(paneth) <- 'gs'

fpkm <- d[rownames(d) %in% paneth$gs,]
pseudoc <- 1
lfpkm <- log2(fpkm+pseudoc)

library(pheatmap)
pheatmap(lfpkm, annotation_col = meta)


ave <- colMeans(lfpkm)
paneth_treat <- data.frame(ave =ave, id = names(ave), stringsAsFactors = FALSE)
paneth_treat_m <- paneth_treat[grepl('WT', paneth_treat$id),]
paneth_treat_m$model <- substr(paneth_treat_m$id, 0, 7)
fc <- function(model, data, col1, col2) {
  d <- data[data$model == model,]
  fc <- d[grepl(col1, d$id), 'ave'] - d[grepl(col2, d$id), 'ave']
  return(fc)
}

panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'CTX', 'NO_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))




panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'NO_EGF', 'WT_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))




panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m, 'CTX', 'WT_EGF')


scores <- data.frame(row.names=names(panethIndScore), PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% c('CRC1139', 'CRC1502', 'CRC1620'),  'ctx'] <- 'R'
#scores$lPIS <- log(scores$PIS)/log(2)
scores$lPIS <- scores$PIS
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
