
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]

r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'



pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))

d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/fpkm.tsv.gz'), sep="\t", header=T)
paneth <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/paneth_s_r_signature'), sep="\t", header=F)
wanted <- c('ATOH1','DEFA5','DEFA6','DLL1','GFI1','AREG', 'EREG','EGF','HBEGF','TGFA','BTC')
h_wanted <- paste0('H_', wanted)

ppaneth <- data.frame(gs=wanted, hgs=h_wanted, class=c(rep('paneth core',5), rep('ligands',6)))


colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)
paneth$class <- 'extended'
paneth <- paneth[!paneth$gs %in% ppaneth$gs,]

paneth <- rbind(paneth, ppaneth)
fpkm <- d[rownames(d) %in% paneth$hgs,]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]
fpkm_pdo_basali <- fpkm[,colnames(fpkm) %in% pdo_basali$id]
pseudoc <- 0.0001
lfpkm_pdo_basali <- log(fpkm_pdo_basali+pseudoc)
lfpkm_pdo_treat <- log(fpkm_pdo_treat+pseudoc)
rownames(lfpkm_pdo_basali) <- substr(rownames(lfpkm_pdo_basali), 3, nchar(rownames(lfpkm_pdo_basali)))
rownames(lfpkm_pdo_treat) <- substr(rownames(lfpkm_pdo_treat), 3, nchar(rownames(lfpkm_pdo_treat)))
pheatmap(lfpkm_pdo_basali, show_rownames=TRUE, show_colnames=FALSE)
pdo_treat_annot <- pdo_treat[,c('id','type','ctx')]
pdo_treat_annot$type <- sapply(strsplit(as.character(pdo_treat_annot$type), '.', fixed=T), function(x){x[1]})
rownames(pdo_treat_annot) <- pdo_treat_annot$id
pdo_treat_annot$id <- NULL
pdo_treat_annot$model <- substr(rownames(pdo_treat_annot), 0, 7)
pdo_treat_annot <- pdo_treat_annot[order(pdo_treat_annot$ctx, pdo_treat_annot$model, pdo_treat_annot$type) ,]
lfpkm_pdo_treat <- lfpkm_pdo_treat[, match(rownames(pdo_treat_annot), colnames(lfpkm_pdo_treat))]
pheatmap(lfpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F)

annot_row <- data.frame(row.names=paneth$gs, class=paneth$class)
pheatmap(lfpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F, annotation_row=annot_row)

lfpkm_pdo_treat <- lfpkm_pdo_treat[ match(rownames(annot_row), rownames(lfpkm_pdo_treat)),]
pheatmap(lfpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F, cluster_rows=F, annotation_row=annot_row)

ave <- colMeans(fpkm_pdo_treat)
paneth_treat <- data.frame(ave =ave, id = names(ave))
paneth_treat_m <- merge(pdo_treat, paneth_treat, by="id")
paneth_treat_m <- paneth_treat_m[order(paneth_treat_m$model, paneth_treat_m$type),]

fc <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[1, 'ave'] / d[2, 'ave'], d[4, 'ave'] / d[3, 'ave'])) # due to order
    } else {
    fc <- d[1, 'ave'] / d[2, 'ave']
  }
  return(fc)
}

paneth <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m)

basal <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, 'ave'], d[3, 'ave'])) # due to order
  } else {
    fc <- d[1, 'ave'] / d[2, 'ave'] # bug here
  }
  return(fc)
}

basald <- sapply(unique(paneth_treat_m$model), basal, paneth_treat_m)
scores <- data.frame(row.names=names(basald), PNS=basald, PIS=paneth)

scores$model <- rownames(scores)
scores$ctx <- 'S'
r <- unique(substr(r$id,0,7))
scores[scores$model %in% r, 'ctx'] <- 'R'
scores$lPIS <- log(scores$PIS)/log(2)
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
ggplot(scores, aes(y=PNS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))

ave <- colMeans(fpkm_pdo_basali)
paneth_basali <- data.frame(ave =ave, id = names(ave))
paneth_basali$model <- substr(paneth_basali$id, 0,7)

basal_basal <- function(model, data) {
  d <- data[data$model == model,]
  return(mean(d$ave))
}

d_basal_basal <- sapply(unique(paneth_basali$model), basal_basal, paneth_basali)
b_scores <- data.frame(row.names=names(d_basal_basal), PNS=d_basal_basal)

b_scores$model <- rownames(b_scores)
ire <- c('CRC0534','CRC1272','CRC0076','CRC0059','CRC0252','CRC0297','CRC0069','CRC0542')
b_scores$ire <- 'no'
b_scores[b_scores$model %in% ire, 'ire'] <- 'yes'
ba <- c('CRC0022','CRC0066','CRC0076','CRC0177','CRC0515','CRC0475')
b_scores$ba <- 'no'; b_scores[b_scores$model %in% ba, 'ba'] <- 'yes'
b_scores$name <- b_scores$model
b_scores[b_scores$ire =="no" & b_scores$ba == "no",'name']<- ''
table(b_scores$name)


ggplot(b_scores, aes(y=PNS,x=reorder(model, -PNS), fill=ire, color=ba))+geom_col()+ylab("basal PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=9, angle = 90, vjust=0.5, hjust=1))+scale_fill_manual(values=c("grey","red"))+scale_color_manual(values=c('grey','black'))+scale_x_discrete(breaks=b_scores$name,labels=b_scores$name)

### heatmap with FC 

mat <- t(fpkm_pdo_treat+pseudoc)
colnames(mat) <- substr(colnames(mat), 3, nchar(colnames(mat)))
mmat <- merge(pdo_treat, mat, by.x="id", by.y="row.names")
mmat <- mmat[order(mmat$model, mmat$type),]

gene_fc <- function(gene, totdata) {
  fcsplit <- function(model, data) {
    d <- data[data$model == model,]
    if (nrow(d) == 4) {
      fc <- c(d[1, gene] / d[2, gene], d[4, gene] / d[3, gene]) # due to order
    } else {
      fc <- c(d[1, gene] / d[2, gene])
    }
    return(fc)
  }
  
  mnames <- function(model, data) {
    d <- data[data$model == model,]
    if (nrow(d) == 4) {
      fc <- c(model, model)
    } else {
      fc <- model
    }
    return(fc)
  }
  
  
  fc <- sapply(unique(totdata$model), fcsplit, totdata)
  fc_unl <- do.call(c, fc)
  names(fc_unl) <- unlist(sapply(unique(totdata$model), mnames, totdata))
  fc_scores <- data.frame(model=names(fc_unl), FC=fc_unl)
  fc_scores
}

all_fc <- lapply(unique(rownames(lfpkm_pdo_treat)), gene_fc, mmat) # no merge on genes to get matrix! TODO

tmp <- all_fc[[1]]
all_all_fc <- do.call(cbind, all_fc)
all_all_fc <- all_all_fc[, !grepl('model', colnames(all_all_fc), fixed=TRUE)]

colnames(all_all_fc) <- unique(rownames(lfpkm_pdo_treat))
rownames(all_all_fc) <- make.unique(as.character(tmp$model))
annot <- data.frame(row.names=rownames(all_all_fc), model=tmp$model)
annot$ctx <- 'S'
annot[annot$model %in% r, 'ctx'] <- 'R'

annot <- annot[order(annot$ctx, annot$model),]
tall_all_fc <- t(all_all_fc)
tall_all_fc <- tall_all_fc[, match(rownames(annot), colnames(tall_all_fc))]

#tall_all_fc[is.infinite(tall_all_fc)]  <- 200
#tall_all_fc[is.nan(tall_all_fc)]  <- 0.000001

#tall_all_fc[tall_all_fc==0]  <- 0.00001

pheatmap(log2(tall_all_fc), show_rownames=TRUE, show_colnames=FALSE, annotation_col = annot, cluster_cols = F)


rm <- rowMeans(tall_all_fc)
tall <- tall_all_fc[order(-rm),]

g <- c('DEFA6','DLL1','ATOH1','GFI1')
annot_row <- data.frame(row.names=rownames(tall))
annot_row$smallsign <- 'no'
annot_row[rownames(annot_row) %in% g, 'smallsign'] <- 'yes'

pheatmap(log2(tall), show_rownames=TRUE, show_colnames=TRUE, annotation_col = annot, annotation_row=annot_row, cluster_cols = F, cluster_rows=F)

data <- data.frame(mean=log2(colMeans(tall)), model=substr(colnames(tall),0,7))
data$ctx <- 'S'
data$x <- colnames(tall)
data[data$model %in% r, 'ctx'] <- 'R'
ggplot(data, aes(y=mean,x=reorder(x, -mean), fill=ctx,))+geom_col()+ylab("average FC after")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=9, angle = 90, vjust=0.5, hjust=1))+scale_fill_manual(values=c("red","blue"))

########
g <- c('DEFA6','DLL1','ATOH1','GFI1')

annot_col <- data.frame(row.names=rownames(lfpkm_pdo_treat))
annot_col$smallsign <- 'no'
annot_col[rownames(annot_col) %in% g, 'smallsign'] <- 'yes'
annot_col <- annot_col[order(annot_col$smallsign),, drop=FALSE]
lfpkm_pdo_treat2 <- lfpkm_pdo_treat[match(rownames(annot_col), rownames(lfpkm_pdo_treat)),]
pheatmap(lfpkm_pdo_treat2, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F, annotation_row=annot_col, cluster_rows=F)


## fpkm
#> dim(d)
#[1] 34466  1225
### tmm
d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/tmm.tsv.gz'), sep="\t", header=T)
dim(d)



meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]

r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'



pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))

paneth <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/scRNA/local/share/data/paneth_s_r_signature'), sep="\t", header=F)

colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)
fpkm <- d[rownames(d) %in% paneth$hgs,]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]
fpkm_pdo_basali <- fpkm[,colnames(fpkm) %in% pdo_basali$id]
pseudoc <- 0.0001
lfpkm_pdo_basali <- log(fpkm_pdo_basali+pseudoc)
lfpkm_pdo_treat <- log(fpkm_pdo_treat+pseudoc)
rownames(lfpkm_pdo_basali) <- substr(rownames(lfpkm_pdo_basali), 3, nchar(rownames(lfpkm_pdo_basali)))
rownames(lfpkm_pdo_treat) <- substr(rownames(lfpkm_pdo_treat), 3, nchar(rownames(lfpkm_pdo_treat)))
pheatmap(lfpkm_pdo_basali, show_rownames=TRUE, show_colnames=FALSE)
pdo_treat_annot <- pdo_treat[,c('id','type','ctx')]
pdo_treat_annot$type <- sapply(strsplit(as.character(pdo_treat_annot$type), '.', fixed=T), function(x){x[1]})
rownames(pdo_treat_annot) <- pdo_treat_annot$id
pdo_treat_annot$id <- NULL
pdo_treat_annot$model <- substr(rownames(pdo_treat_annot), 0, 7)
pdo_treat_annot <- pdo_treat_annot[order(pdo_treat_annot$ctx, pdo_treat_annot$model, pdo_treat_annot$type) ,]
lfpkm_pdo_treat <- lfpkm_pdo_treat[, match(rownames(pdo_treat_annot), colnames(lfpkm_pdo_treat))]
pheatmap(lfpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F)

ave <- colMeans(fpkm_pdo_treat)
paneth_treat <- data.frame(ave =ave, id = names(ave))
paneth_treat_m <- merge(pdo_treat, paneth_treat, by="id")
paneth_treat_m <- paneth_treat_m[order(paneth_treat_m$model, paneth_treat_m$type),]

fc <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[1, 'ave'] / d[2, 'ave'], d[4, 'ave'] / d[3, 'ave'])) # due to order
  } else {
    fc <- d[1, 'ave'] / d[2, 'ave']
  }
  return(fc)
}

paneth <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m)

basal <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, 'ave'], d[3, 'ave'])) # due to order
  } else {
    fc <- d[1, 'ave'] / d[2, 'ave']
  }
  return(fc)
}

basald <- sapply(unique(paneth_treat_m$model), basal, paneth_treat_m)
scores <- data.frame(row.names=names(basald), PNS=basald, PIS=paneth)

scores$model <- rownames(scores)
scores$ctx <- 'S'
r <- unique(substr(r$id,0,7))
scores[scores$model %in% r, 'ctx'] <- 'R'
scores$lPIS <- log(scores$PIS)/log(2)
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
ggplot(scores, aes(y=PNS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))

ave <- colMeans(fpkm_pdo_basali)
paneth_basali <- data.frame(ave =ave, id = names(ave))
paneth_basali$model <- substr(paneth_basali$id, 0,7)

basal_basal <- function(model, data) {
  d <- data[data$model == model,]
  return(mean(d$ave))
}

d_basal_basal <- sapply(unique(paneth_basali$model), basal_basal, paneth_basali)
b_scores <- data.frame(row.names=names(d_basal_basal), PNS=d_basal_basal)

b_scores$model <- rownames(b_scores)
ire <- c('CRC0534','CRC1272','CRC0076','CRC0059','CRC0252','CRC0297','CRC0069','CRC0542')
b_scores$ire <- 'no'
b_scores[b_scores$model %in% ire, 'ire'] <- 'yes'
ba <- c('CRC0022','CRC0066','CRC0076','CRC0177','CRC0515','CRC0475')
b_scores$ba <- 'no'; b_scores[b_scores$model %in% ba, 'ba'] <- 'yes'
b_scores$name <- b_scores$model
b_scores[b_scores$ire =="no" & b_scores$ba == "no",'name']<- ''
table(b_scores$name)


ggplot(b_scores, aes(y=PNS,x=reorder(model, -PNS), fill=ire, color=ba))+geom_col()+ylab("basal PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=9, angle = 90, vjust=0.5, hjust=1))+scale_fill_manual(values=c("grey","red"))+scale_color_manual(values=c('grey','black'))+scale_x_discrete(breaks=b_scores$name,labels=b_scores$name)


### heatmap with FC 

mat <- t(fpkm_pdo_treat+pseudoc)
colnames(mat) <- substr(colnames(mat), 3, nchar(colnames(mat)))
mmat <- merge(pdo_treat, mat, by.x="id", by.y="row.names")
mmat <- mmat[order(mmat$model, mmat$type),]

gene_fc <- function(gene, totdata) {
  fcsplit <- function(model, data) {
    d <- data[data$model == model,]
    if (nrow(d) == 4) {
      fc <- c(d[1, gene] / d[2, gene], d[4, gene] / d[3, gene]) # due to order
    } else {
      fc <- c(d[1, gene] / d[2, gene])
    }
    return(fc)
  }
  
  mnames <- function(model, data) {
    d <- data[data$model == model,]
    if (nrow(d) == 4) {
      fc <- c(model, model)
    } else {
      fc <- model
    }
    return(fc)
  }
  
  
  fc <- sapply(unique(totdata$model), fcsplit, totdata)
  fc_unl <- do.call(c, fc)
  names(fc_unl) <- unlist(sapply(unique(totdata$model), mnames, totdata))
  fc_scores <- data.frame(model=names(fc_unl), FC=fc_unl)
  fc_scores
}

all_fc <- lapply(unique(rownames(lfpkm_pdo_treat)), gene_fc, mmat) # no merge on genes to get matrix! TODO

tmp <- all_fc[[1]]
all_all_fc <- do.call(cbind, all_fc)
all_all_fc <- all_all_fc[, !grepl('model', colnames(all_all_fc), fixed=TRUE)]

colnames(all_all_fc) <- unique(rownames(lfpkm_pdo_treat))
rownames(all_all_fc) <- make.unique(as.character(tmp$model))
annot <- data.frame(row.names=rownames(all_all_fc), model=tmp$model)
annot$ctx <- 'S'
annot[annot$model %in% r, 'ctx'] <- 'R'

annot <- annot[order(annot$ctx, annot$model),]
tall_all_fc <- t(all_all_fc)
tall_all_fc <- tall_all_fc[, match(rownames(annot), colnames(tall_all_fc))]

#tall_all_fc[is.infinite(tall_all_fc)]  <- 200
#tall_all_fc[is.nan(tall_all_fc)]  <- 0.000001

#tall_all_fc[tall_all_fc==0]  <- 0.00001

pheatmap(log2(tall_all_fc), show_rownames=TRUE, show_colnames=FALSE, annotation_col = annot, cluster_cols = F)


rm <- rowMeans(tall_all_fc)
tall <- tall_all_fc[order(-rm),]

g <- c('DEFA6','DLL1','ATOH1','GFI1')
annot_row <- data.frame(row.names=rownames(tall))
annot_row$smallsign <- 'no'
annot_row[rownames(annot_row) %in% g, 'smallsign'] <- 'yes'

pheatmap(log2(tall), show_rownames=TRUE, show_colnames=TRUE, annotation_col = annot, annotation_row=annot_row, cluster_cols = F, cluster_rows=F)

data <- data.frame(mean=log2(colMeans(tall)), model=substr(colnames(tall),0,7))
data$ctx <- 'S'
data$x <- colnames(tall)
data[data$model %in% r, 'ctx'] <- 'R'
ggplot(data, aes(y=mean,x=reorder(x, -mean), fill=ctx,))+geom_col()+ylab("average FC after")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=9, angle = 90, vjust=0.5, hjust=1))+scale_fill_manual(values=c("red","blue"))



########### indagini IRe
wanted <- c('ATOH1','DEFA5','DEFA6','DLL1','GFI1')
h_wanted <- paste0('H_', wanted)

all_all_fc[grepl('CRC0069', rownames(all_all_fc)), colnames(all_all_fc) %in% wanted]
all_all_fc[grepl('CRC0322', rownames(all_all_fc)), colnames(all_all_fc) %in% wanted]
