meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]

r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'

pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))

d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)
#wanted <- data.frame(gs=c('ATOH1','DEFA5','DEFA6','DLL1','GFI1','AREG', 'EREG','EGF','HBEGF','TGFA','BTC'))
paneth <- data.frame(gs=c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4"))

colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)

fpkm <- d[rownames(d) %in% paneth$hgs,]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]
fpkm_pdo_basali <- fpkm[,colnames(fpkm) %in% pdo_basali$id] #no loss of replicates: all . here
pseudoc <- 1
lfpkm_pdo_basali <- log(fpkm_pdo_basali+pseudoc)
lfpkm_pdo_treat <- log(fpkm_pdo_treat+pseudoc)

ave <- colMeans(lfpkm_pdo_treat) # run this way, changes a lot
#ave <- colMeans(fpkm_pdo_treat) 
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

panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m)


basal <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, 'ave'], d[3, 'ave'])) # due to order
  } else {
    fc <- d[2, 'ave']
  }
  return(fc)
}

basald <- sapply(unique(paneth_treat_m$model), basal, paneth_treat_m)
scores <- data.frame(row.names=names(basald), PNS=basald, PIS=panethIndScore)

scores$model <- rownames(scores)
scores$ctx <- 'S'
r <- unique(substr(r$id,0,7))
scores[scores$model %in% r, 'ctx'] <- 'R'
scores$lPIS <- log(scores$PIS)/log(2)
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
#ggsave('lPIS_tmm.svg')
ggplot(scores, aes(y=PNS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
#ggsave('PNS_tmm.svg')
#write.table(scores, file='PIS_scores_tmm.tsv', sep="\t", quote=FALSE, row.names=F)

pdo <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/cetuximab/atp/merge.tsv', sep='\t', header=T)
mscores <- merge(scores, pdo, by.y="case", by.x='model')

write.table(mscores, file='PIS_scores_atp_tmm.tsv', sep="\t", quote=FALSE, row.names=F)
ggplot(mscores, aes(y=lPIS,x=reorder(model, -lPIS),fill=CTG_5000))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_distiller(palette="RdYlBu")
ggsave('PNS_tmm_gradientviability.svg')

#### PIS vs ETS1 basale

efpkm <- d[rownames(d) == "H_ETS1",]
efpkm_pdo_treat <- efpkm[,colnames(efpkm) %in% pdo_treat$id]
efpkm_pdo_basali <- efpkm[,colnames(efpkm) %in% pdo_basali$id] #no loss of replicates: all . here
pseudoc <- 1
elfpkm_pdo_basali <- log(efpkm_pdo_basali+pseudoc)
elfpkm_pdo_treat <- log(efpkm_pdo_treat+pseudoc)
ets <- data.frame(id=names(elfpkm_pdo_treat), ETS1=as.numeric(elfpkm_pdo_treat))
ets_m <- merge(pdo_treat, ets, by="id")
ets_m <- ets_m[order(ets_m$type),]


efc <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    #fc <- mean(c(d[1, 'ave'] / d[2, 'ave'], d[4, 'ave'] / d[3, 'ave'])) # due to order
    res <- mean(d[2,'ETS1'],d[3,'ETS1'])
  } else {
    #fc <- d[1, 'ave'] / d[2, 'ave']
    res <- d[2,'ETS1']
  }
  return(res)
}

ets_nt <- as.data.frame(sapply(unique(ets_m$model), efc, ets_m))

mets_pis <- merge(ets_nt, mscores, by.x="row.names", by.y="model")
ggplot(mets_pis, aes(x=lPIS, y=x, color=CTG_5000))+geom_point()+geom_smooth(method="lm")+theme_bw()+ylab('ETS1 noEGF')+scale_color_distiller(palette="RdYlBu")

etsEGF <- data.frame(id=names(elfpkm_pdo_basali), ETS1=as.numeric(elfpkm_pdo_basali))
etsEGF$model <- substr(etsEGF$id, 0,7)

egf_ave <- function(model, data) {
  d <- data[data$model == model,]
  return(mean(d$ETS1))
}

ets_basale <- as.data.frame(sapply(unique(etsEGF$model), egf_ave, etsEGF))
colnames(ets_basale) <- 'ETS'
etsb_pis <- merge(ets_basale, mscores, by.x="row.names", by.y="model")
ggplot(etsb_pis, aes(x=PIS, y=ETS, color=CTG_5000))+geom_point()+geom_smooth(method="lm")+theme_bw()+ylab('ETS1 EGF')+scale_color_distiller(palette="RdYlBu")


ggplot(mets_pis, aes(x=PNS, y=x, color=CTG_5000))+geom_point()+geom_smooth(method="lm")+theme_bw()+ylab('ETS1 noEGF')+scale_color_distiller(palette="RdYlBu")

##########

pdo_treat_annot <- pdo_treat[,c('id','type','ctx')]
pdo_treat_annot$type <- sapply(strsplit(as.character(pdo_treat_annot$type), '.', fixed=T), function(x){x[1]})
rownames(pdo_treat_annot) <- pdo_treat_annot$id
pdo_treat_annot$id <- NULL
pdo_treat_annot$model <- substr(rownames(pdo_treat_annot), 0, 7)
pdo_treat_annot <- pdo_treat_annot[order(pdo_treat_annot$ctx, pdo_treat_annot$model, pdo_treat_annot$type) ,]
lfpkm_pdo_treat <- lfpkm_pdo_treat[, match(rownames(pdo_treat_annot), colnames(lfpkm_pdo_treat))]
pheatmap(lfpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F)
pheatmap(fpkm_pdo_treat, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F)


###  metagene at the fc level
single_fc <- function(model, data, info) {
  d <- t(data[,grepl(model, colnames(data))])
  m <- merge(info, d, by.y="row.names", by.x="id")
  m <- m[order(m$type),]
  expr <- m[, seq(8, ncol(m))]
  expr[expr==0] <- 0.00001
  if (nrow(expr) == 4) {
    fc <- (expr[1,] / expr[2,] +  expr[4,] / expr[3,])/2 # due to order
  } else {
    fc <- expr[1,] / expr[2,]
  }
  return(fc)
}

singlegene_indScore <- lapply(unique(paneth_treat_m$model), single_fc, fpkm_pdo_treat, pdo_treat)
sis <- do.call(rbind,singlegene_indScore)
rownames(sis) <- unique(paneth_treat_m$model)
sscores <- rowMeans(log(sis)/log(2))

scores <- data.frame(row.names=names(sscores), PIS=sscores)

scores$model <- rownames(scores)
scores$ctx <- 'S'
scores[scores$model %in% r, 'ctx'] <- 'R'
ggplot(scores, aes(y=lPIS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))
#ggsave('lPIS_tmm.svg')
ggplot(scores, aes(y=PNS,x=reorder(model, -lPIS),fill=ctx))+geom_col()+ylab("PNS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))


ggplot(scores, aes(y=PIS,x=reorder(model, -PIS),fill=ctx))+geom_col()+ylab("PIS")+xlab("Model")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("red","blue"))


##########3 correlations between samples ############
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]
r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'
pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))
tmm <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)

tmm_pdo_treat <- tmm[,colnames(tmm) %in% pdo_treat$id]
tmm_pdo_basali <- tmm[,colnames(tmm) %in% pdo_basali$id] #no loss of replicates: all . here
pseudoc <- 1
ltmm_pdo_basali <- log(tmm_pdo_basali+pseudoc)
ltmm_pdo_treat <- log(tmm_pdo_treat+pseudoc)

sds <- apply(ltmm_pdo_treat, 1, sd)
sds <- sds[order(-sds)]
topvar <- ltmm_pdo_treat[rownames(ltmm_pdo_treat) %in% names(head(sds, n=3000)),]

cors <- cor(topvar)
annot <- data.frame(row.names=pdo_treat$id, model=pdo_treat$model, type=substr(pdo_treat$type,0,17), ctx=pdo_treat$ctx)
ann_colors = list(
  ctx = c(S="blue", R="firebrick")
  #CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  #GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)
pheatmap(cors, annotation_row = annot, show_rownames = FALSE, show_colnames=TRUE, annotation_colors=ann_colors)
cors[grepl('CRC0322',rownames(cors)),grepl('CRC0322',colnames(cors))]
upper <- cors[upper.tri(cors)]
pd <- data.frame(pearson=upper)

v <- 0.9824710  
ggplot(data=pd, aes(pearson))+geom_histogram(bins=30, fill="white",color="black")+theme_bw()+theme(axis.text=element_text(size=20))+geom_vline(xintercept=v, color="red",size=1)
coors <- sapply(models, function(x) {y <- cors[grepl(x,rownames(cors)),grepl(x,colnames(cors))]; y[upper.tri(y)]})

pd <- data.frame(pearson=unlist(coors))
ggplot(data=pd, aes(pearson))+geom_histogram(bins=30, fill="white",color="black")+theme_bw()+theme(axis.text=element_text(size=20))+geom_vline(xintercept=v, color="red",size=1)
