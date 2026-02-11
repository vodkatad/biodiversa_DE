library(ggplot2)
setwd('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/chemio_collection')
#d <- read.table('tmm.tsv.gz', sep="\t", header=T)
d0 <- read.table('vsd.tsv.gz', sep="\t", header=T)

s <- read.table('old_samples.tsv', sep="\t", header=F)

d <- d0[grepl('M_', rownames(d0)),]
d <- d[, s$V2]

dm <- d0[grepl('H_', rownames(d0)),]
dm <- dm[, s$V2]

#rm <- rowMeans(d)
#d <- d[rm > 5,]
alle <- unlist(d)
pd <- data.frame(expr=alle, row.names=names(alle))
pd$type <- 'M'

#rm <- rowMeans(dm)
#dm <- dm[rm > 5,]
allem <- unlist(dm)
pdm <- data.frame(expr=allem, row.names=names(allem))
pdm$type <- 'H'
pd <- rbind(pd, pdm)
#ggplot(data=pd, aes(x=log2(expr+1)))+geom_histogram(bins=50)
ggplot(data=pd, aes(x=expr, fill=type))+geom_histogram(bins=50, position='dodge')
# vsd <- pd

cl <- read.table('/mnt/cold1/bioinfotree/task/gencode/dataset/mmusculus_hsapiens_combined/M16_27/gene_class.tsv', sep="\t", header=F)

cl_s <- cl[cl$V2 %in% rownames(dm),]
cl_s <- cl_s[match(rownames(dm), cl_s$V2),]

pc <- unlist(dm[cl_s$V1=="protein_coding",])
mi <- min(pc)
pcd <- data.frame(expr=pc)
pcd$type <- 'pc'
#ggplot(data=pcd, aes(x=log2(expr+1)))+geom_histogram(bins=50)
ggplot(data=pcd, aes(x=expr))+geom_histogram(bins=50)

pc1 <- unlist(d[cl_s$V1!="protein_coding",])
pcd1 <- data.frame(expr=pc1)
pcd1$type <- 'npc'
pcd <- rbind(pcd, pcd1)
ggplot(data=pcd, aes(x=expr, fill=type))+geom_histogram(bins=50, position="dodge")
