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