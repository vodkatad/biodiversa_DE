
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]

r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'



pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))

fpkm <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm.tsv.gz'), sep="\t", header=T)

fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]

pdo_treat_annot <- pdo_treat[,c('id','type','ctx')]
pdo_treat_annot$type <- sapply(strsplit(as.character(pdo_treat_annot$type), '.', fixed=T), function(x){x[1]})
rownames(pdo_treat_annot) <- pdo_treat_annot$id
pdo_treat_annot$id <- NULL
pdo_treat_annot$model <- substr(rownames(pdo_treat_annot), 0, 7)
pdo_treat_annot <- pdo_treat_annot[order(pdo_treat_annot$ctx, pdo_treat_annot$model, pdo_treat_annot$type) ,]


wanted <- c('ATOH1','DEFA5','DEFA6','DLL1','GFI1','AREG', 'EREG','EGF','HBEGF','TGFA','BTC', 'ETS1')
h_wanted <- paste0('H_', wanted)
ppaneth <- data.frame(gs=wanted, hgs=h_wanted, class=c(rep('paneth core',5), rep('ligands',6), 'new friend'), stringsAsFactors = F)

fpkm_paneth <- fpkm[rownames(fpkm) %in% ppaneth$hgs,]
fpkm_pdo_treat_paneth <- fpkm_paneth[,colnames(fpkm_paneth) %in% pdo_treat$id]
pseudoc <- 1
lfpkm_pdo_treat_paneth <- log(fpkm_pdo_treat_paneth+pseudoc)
lfpkm_pdo_treat_paneth <- lfpkm_pdo_treat_paneth[, match(rownames(pdo_treat_annot), colnames(lfpkm_pdo_treat_paneth))]
pheatmap(lfpkm_pdo_treat_paneth, show_rownames=TRUE, show_colnames=FALSE, annotation_col = pdo_treat_annot, cluster_cols = F)

cols <- apply(sapply(c('CRC0327', 'CRC0542', 'CRC0322'), function(x) { grepl(x, colnames(lfpkm_pdo_treat_paneth))}), 1, any)
keep <- t(lfpkm_pdo_treat_paneth[rownames(lfpkm_pdo_treat_paneth)=="H_ETS1", cols])

mkeep <- merge(keep, pdo_treat_annot, by="row.names")
mkeep$x <- seq(1, nrow(mkeep))
ggplot(data=mkeep, aes(x=x, fill=model, alpha=type, y=H_ETS1))+geom_col(position="dodge")+scale_alpha_discrete(range = c(0.8, 0.4))


tfpkm <- t(fpkm_pdo_treat)
mkeep <- merge(pdo_treat_annot, tfpkm, by="row.names")

write.table(mkeep, file="~/pdo_treat_matrix.tsv", sep="\t", quote=FALSE)
