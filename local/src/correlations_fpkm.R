setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/batch5_h')
fpkm_h <- read.table(gzfile('fpkm.tsv.gz'), sep="\t", header=T)
fpkm_mh <- read.table(gzfile('../batch5_hm/fpkm.tsv.gz'), sep="\t", header=T)
library(tidyverse)
rownames(fpkm_mh) <- str_remove(rownames(fpkm_mh),"H_")
common <- intersect(rownames(fpkm_h), rownames(fpkm_mh))
ff_h <- fpkm_h[rownames(fpkm_h) %in% common,]
ff_mh <- fpkm_mh[rownames(fpkm_mh) %in% common,]
ff_h <- ff_h[order(rownames(ff_h)),]
ff_mh <- ff_mh[order(rownames(ff_mh)),]
all(rownames(ff_h)==rownames(ff_mh))
all(colnames(ff_h)==colnames(ff_mh))

ti1<- Sys.time(); correlations <- cor(t(ff_mh), t(ff_h)); ti2 <- Sys.time(); print(ti2-ti1)
# 7'
cors <- diag(correlations)
histogram(cors, xlim=c(0.9,1), breaks=seq(0.9,1,by=0.01))
---
fpkm_h <- read.table(gzfile('counts.tsv.gz'), sep="\t", header=T, row.names=1)
fpkm_mh <- read.table(gzfile('../batch5_hm/counts.tsv.gz'), sep="\t", header=T, row.names=1)
rownames(fpkm_mh) <- str_remove(rownames(fpkm_mh),"H_")
common <- intersect(rownames(fpkm_h), rownames(fpkm_mh))
ff_h <- fpkm_h[rownames(fpkm_h) %in% common,]
ff_mh <- fpkm_mh[rownames(fpkm_mh) %in% common,]
ff_h <- ff_h[order(rownames(ff_h)),]
ff_mh <- ff_mh[order(rownames(ff_mh)),]
all(rownames(ff_h)==rownames(ff_mh))
all(colnames(ff_h)==colnames(ff_mh))

ti1<- Sys.time(); ccorrelations <- cor(t(ff_mh), t(ff_h)); ti2 <- Sys.time(); print(ti2-ti1)
# 7'
ccors <- diag(ccorrelations)
summary(ccors)
histogram(ccors, xlim=c(0.9,1), breaks=seq(0.9,1,by=0.01))

  
  
