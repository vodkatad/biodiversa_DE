library(ggplot2)
expr_f <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd_H.tsv.gz'
ss_f <-  '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data'
meta <- c('ATOH1', 'DLL1', 'GFI1', 'DEFA5', 'DEFA6')
single <- 'HES1'
# metagene Paneth basale e HES1 vs risposta cetux PDX

expr <- read.table(gzfile(expr_f), sep="\t", header=T, stringsAsFactors = F)
ss <- read.table(ss_f, sep="\t", header=T, stringsAsFactors = F)

wanted <- 'LMX_BASALE'
wanted_lmodel <- ss[grepl(wanted, ss$type), 'id']
expr_w <- expr[, colnames(expr) %in% wanted_lmodel]

sm <- substr(colnames(expr_w), 0, 7)

smodel <- unique(sm)

ave_models <- function(model, my_expr) {
  my_sel_expr <- my_expr[, grepl(model, colnames(my_expr)), drop=FALSE]
  if (ncol(my_sel_expr) > 1) {
    return(as.data.frame(rowMeans(my_sel_expr)))
  } else {
    return(my_sel_expr)
  }
}

ave_expr <- lapply(smodel, ave_models, expr_w)
ave_expr_df <- as.data.frame(do.call(cbind, ave_expr))
rownames(ave_expr_df) <- rownames(expr_w)
colnames(ave_expr_df) <- smodel

cetuxi <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/last_march2024/last_cet_march2024.txt', sep="\t", header=TRUE, stringsAsFactors = F)
colnames(cetuxi)[2] <- 'Cetuxi_VolVar3WKS'
colnames(cetuxi)[3] <- 'Cetuxi_VolVar6WKS'
cetuxi[2] <- cetuxi[2]*100
cetuxi[3] <- cetuxi[3]*100


meta_ave_expr_df <- ave_expr_df[rownames(ave_expr_df) %in% meta, ]
meta_df <- colMeans(meta_ave_expr_df)
single_ave_expr_df <- ave_expr_df[rownames(ave_expr_df) == single, , drop=F]

## ssGSEA scores for paneth?

data <- t(rbind(meta_df, single_ave_expr_df))
colnames(data)[1] <- 'meta_paneth'

#sdata <- scale(data)
sdata <- data
m <- merge(sdata, cetuxi, by.x='row.names', by.y='CASE')
ggplot(data=m, aes(x=HES1, y=Cetuxi_VolVar3WKS))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size=15)
ggplot(data=m, aes(x=meta_paneth, y=Cetuxi_VolVar3WKS))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size=15)

cor.test(m$HES1, m$Cetuxi_VolVar3WKS)
cor.test(m$meta_paneth, m$Cetuxi_VolVar3WKS)


## ssGSEA scores on basali only?

## PIS revival ###############################################
meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", header=T)
pdo_basali <- meta[grepl('LMO_BASALE', meta$type),]
pdo_treat <- meta[grepl('LMO_cetuxi', meta$type),]

r <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', header=T)
pdo_treat$ctx <- 'S'
pdo_treat[pdo_treat$id %in% r$id, 'ctx'] <- 'R'

pdo_treat$model <- substr(pdo_treat$id, 0,7)
basali <- unique(substr(pdo_basali[,'id'],0,7))

d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)
paneth <- data.frame(gs=c('ATOH1','DEFA5','DEFA6','DLL1','GFI1'))

colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)

fpkm <- d[rownames(d) %in% paneth$hgs,]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]
fpkm_pdo_basali <- fpkm[,colnames(fpkm) %in% pdo_basali$id] #no loss of replicates: all . here
pseudoc <- 1
lfpkm_pdo_basali <- log2(fpkm_pdo_basali+pseudoc) #changed to log2 the 6/12/2023 to be in line with Ire's thesis, results are ==
lfpkm_pdo_treat <- log2(fpkm_pdo_treat+pseudoc)

ave <- colMeans(lfpkm_pdo_treat) # run this way, changes a lot
#ave <- colMeans(fpkm_pdo_treat) 
paneth_treat <- data.frame(ave =ave, id = names(ave))
paneth_treat_m <- merge(pdo_treat, paneth_treat, by="id")
paneth_treat_m <- paneth_treat_m[order(paneth_treat_m$model, paneth_treat_m$type),]

fc <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[1, 'ave'] - d[2, 'ave'], d[4, 'ave'] - d[3, 'ave'])) # due to order XXX -
  } else {
    fc <- d[1, 'ave'] - d[2, 'ave']
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

ggplot(data=scores, aes(x=ctx, y=PIS))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw(base_size=15)
ggplot(data=scores, aes(x=ctx, y=PNS))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw(base_size=15)

res <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/genetic_cet_res_update0625.tsv', sep="\t", header=T, stringsAsFactors = F)
# 16 sono tutti i modelli PDO con trattamenti e dati RNAseq storici biobanca
scores$fourWT <- ifelse(rownames(scores) %in% res$CASE, 'no', 'yes')



### HES1
fpkm <- d[rownames(d) == "H_HES1",, drop=F]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% pdo_treat$id]
fpkm_pdo_basali <- fpkm[,colnames(fpkm) %in% pdo_basali$id] #no loss of replicates: all . here
pseudoc <- 1
lfpkm_pdo_basali <- log2(fpkm_pdo_basali+pseudoc) #changed to log2 the 6/12/2023 to be in line with Ire's thesis, results are ==
lfpkm_pdo_treat <- log2(fpkm_pdo_treat+pseudoc)

#ave <- colMeans(lfpkm_pdo_treat) # run this way, changes a lot
ave <- lfpkm_pdo_treat
paneth_treat <- data.frame(ave =unlist(ave), id = names(ave))
paneth_treat_m <- merge(pdo_treat, paneth_treat, by="id")
paneth_treat_m <- paneth_treat_m[order(paneth_treat_m$model, paneth_treat_m$type),]



panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m)

basald <- sapply(unique(paneth_treat_m$model), basal, paneth_treat_m)
hes_scores <- data.frame(row.names=names(basald), PNS=basald, PIS=panethIndScore)
hes_scores[order(hes_scores$PIS),] # dove PIS è HES1 induction

colnames(hes_scores) <- c('HES1_NT', 'HES1_ind')

m <- merge(hes_scores, scores, by="row.names")

ggplot(data=m, aes(x=HES1_NT, y=PIS, color=ctx, shape=fourWT))+geom_point()+theme_bw(base_family = 11)
ggplot(data=m, aes(x=HES1_ind, y=PIS, color=ctx, shape=fourWT))+geom_point()+theme_bw(base_family = 11)
cor.test(m$HES1_ind, m$PIS)

ms <- m[m$ctx=="S",]
cor.test(ms$HES1_ind, ms$PIS)
## CRC0069 upregola HES1
