library(ggplot2)
library(pheatmap)
texpr <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)
rownames(texpr) <- substr(rownames(texpr), 3, nchar(rownames(texpr)))
fexpr <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm.tsv.gz'), sep="\t", header=T)
rownames(fexpr) <- substr(rownames(fexpr), 3, nchar(rownames(fexpr)))
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=T)
egf <- read.table("~/EGF_role.txt", sep="\t", header=T)

basali <- meta[grepl('LMO_BASALE', meta$type), ]

tbasali <- texpr[, colnames(texpr) %in%  basali$id]
fbasali <- fexpr[, colnames(fexpr) %in%  basali$id]

ligandi <- c('AREG','EREG','EGF','HBEGF','TGFA','BTC')

tbasali_lig <- tbasali[rownames(tbasali)%in% ligandi,]
fbasali_lig <- fbasali[rownames(fbasali)%in% ligandi,]
modello <- substr(colnames(tbasali_lig),0,7)
all( colnames(tbasali_lig) == colnames(fbasali_lig))
tbasali_lig_pcr <- tbasali_lig[, modello %in% egf$model]
fbasali_lig_pcr <- fbasali_lig[, modello %in% egf$model]

PC <- 1
rownames(egf) <- egf$model
egf <- egf[order(egf$EGFrole),]
egf$model <- as.character(egf$model)
m <- data.frame(m=substr(colnames(tbasali_lig_pcr),0,7))
me <- merge(m, egf, by.y="row.names", by.x ="m", all.x=T)
me <- me[order(me$EGFrole),]
me$model <- make.unique(me$model)
rownames(me) <- me$model
d <- log2(fbasali_lig_pcr+PC)
colnames(d) <- make.unique(substr(colnames(d),0,7))
d <- d[, match(me$model, colnames(d))]
me2 <- me
me2$model <- NULL
pheatmap(as.matrix(d), annotation_col=me2, cluster_cols=F)

d <- log2(tbasali_lig_pcr+PC)
colnames(d) <- make.unique(substr(colnames(d),0,7))
d <- d[, match(me$model, colnames(d))]
me$model <- NULL
pheatmap(as.matrix(d), annotation_col=me2, cluster_cols=F)


tt <- as.data.frame(t(log2(tbasali_lig_pcr)+1))
tt$model <-  substr(rownames(tt),0,7)

tegf <- sapply(unique(tt$model), function(x) { mean(tt[tt$model==x, colnames(tt)=="EGF"])})
dd <- data.frame(rownames=names(tegf), EGF=tegf)
m <- merge(dd, egf, by.x="rownames", by.y="model")
plot(m$EGF, m$EGFrole)

tegf <- sapply(unique(tt$model), function(x) { mean(tt[tt$model==x, colnames(tt)=="EREG"])})
dd <- data.frame(rownames=names(tegf), EREG=tegf)
m <- merge(dd, egf, by.x="rownames", by.y="model")
plot(m$EREG, m$EGFrole)

tegf <- sapply(unique(tt$model), function(x) { mean(tt[tt$model==x, colnames(tt)=="TGFA"])})
dd <- data.frame(rownames=names(tegf), TGFA=tegf)
m <- merge(dd, egf, by.x="rownames", by.y="model")
plot(m$TGFA, m$EGFrole)

tegf <- sapply(unique(tt$model), function(x) { mean(tt[tt$model==x, colnames(tt)=="HBEGF"])})
dd <- data.frame(rownames=names(tegf), HBEGF=tegf)
m <- merge(dd, egf, by.x="rownames", by.y="model")
plot(m$HBEGF, m$EGFrole)