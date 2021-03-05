library(ggplot2)
library(pheatmap)
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S')

deg <- read.table('treat_cutoff0.05-cetuxi.vs.NT.deseq2_sign_genedesc.tsv', comment.char="", quote='', header=T, sep="\t")
rownames(deg) <- deg$gene
deg$gene <- NULL

## heatmap of DEG ########
expr <- read.table(gzfile('fpkm.tsv.gz'), sep="\t", header=T)
rownames(expr) <- substr(rownames(expr), 3, nchar(rownames(expr)))
samples_data <- read.table('samples_data', sep="\t", header=T)
rownames(samples_data) <- samples_data$id
samples_data$id <- NULL

PC <- 0.0001

#### merge
deg_up <- deg[deg$log2FoldChange < 0,]
deg_up$description <- as.character(deg_up$description)
pdo_r <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")
pdx_r <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_R/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")
pdx_s <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")
basali <- read.table('../Biodiversa_up5starOK_basal_NT_PDO/class_cutoff0.05-NT.vs.BASALE.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")

not_paired <- read.table('./batch/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")


prepare_df <- function(df, keep=c('log2FoldChange','padj')) {
  rownames(df) <- substr(rownames(df), 3, nchar(rownames(df)))
  df[, keep]
}


pdo_r <- prepare_df(pdo_r, c('log2FoldChange','pvalue'))
colnames(pdo_r) <- c('log2FC_PDO_R', 'pval_PDO_R')
res <- merge(deg_up, pdo_r, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL

pdx_s <- prepare_df(pdx_s)
colnames(pdx_s) <- c('log2FC_PDX_S', 'padj_PDX_S')
res <- merge(res, pdx_s, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL

pdx_r <- prepare_df(pdx_r)
colnames(pdx_r) <- c('log2FC_PDX_R', 'padj_PDX_R')
res <- merge(res, pdx_r, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL

basali <- prepare_df(basali)
colnames(basali) <- c('log2FC_basal_NT', 'padj_basal_NT')
res <- merge(res, basali, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL


not_paired <- prepare_df(not_paired)
colnames(not_paired) <- c('log2FC_notpaired', 'padj_notpaired')
res <- merge(res, not_paired, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL


##
res <- res[order(-res$log2FoldChange, -res$baseMean,res$padj),]
#write.table(res, gzfile('neg_big_merge.tsv.gz'), sep="\t", quote=F, row.names=T)


############## imbuti
mediane <- median(res$baseMean)
res$filter_PDX_S <- ifelse(!is.na(res$padj_PDX_S) & res$padj_PDX_S < 0.05 & res$log2FC_PDX_S < 0, TRUE, FALSE)
res$filter_PDX_R <- ifelse(!is.na(res$padj_PDX_R) & res$padj_PDX_R >= 0.05, TRUE, FALSE) 
res$filter_expr <- ifelse(res$baseMean > mediane, TRUE, FALSE)
res$allfilter <- apply(res[,grepl('filter_', colnames(res))], 1, all)

buckets <- read.table("../../local/share/data/Supplementary_Table_9_list_bucket.txt", sep="\t", header=T)
buck <- aggregate(buckets$Cancer.Tissue.Type, buckets[, c(2,3,4)], FUN=function(x) {as.character(paste0(unique(x), collapse=","))} , simplify=T)
colnames(buck)[4] <- 'cancer_type'

mres <- merge(res, buck, by.x="row.names", by.y="Target", all.x=TRUE)
rownames(mres) <- mres$Row.names
mres$Row.names <- NULL
write.table(res, gzfile('negdeg_filtersteps_tractability_group.tsv.gz'), sep="\t", quote=F, row.names=T)
selection <- mres[!is.na(mres$Tractability.Group) & mres$Tractability.Group <=2 & mres$allfilter,]                      
write.table(selection, gzfile('negdeg_filtersteps_tractability_group_up2.tsv.gz'), sep="\t", quote=F, row.names=T)


deg_expr <- expr[rownames(expr) %in% rownames(selection),]
samples_data <- samples_data[order(samples_data$sample, samples_data$treat),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
annot_row <- data.frame(row.names=rownames(deg), expr=log(deg$baseMean), log2FoldChange=deg$log2FoldChange)
annot_row <- annot_row[rownames(annot_row) %in% rownames(deg_expr),]
pheatmap(log(deg_expr+PC), show_rownames=TRUE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)


gene_fc <- function(sample, expr, samples_data) {
  meta <- samples_data[samples_data$sample == sample,]
  expr <- expr[,colnames(expr) %in% rownames(meta)]
  ctx <- meta[meta$treat =="cetuxi",]
  nt <- meta[meta$treat =="NT",]
  ave_ctx <- rowMeans(expr[,colnames(expr) %in% rownames(ctx)])
  ave_nt <- rowMeans(expr[,colnames(expr) %in% rownames(nt)])
  return(ave_ctx/ave_nt)
}

us <- unique(samples_data$sample)

samplefc <- sapply(us, gene_fc, deg_expr, samples_data)
colnames(samplefc) <- us
pheatmap(samplefc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F)

lfc <- log2(samplefc)
minv <- min(lfc)
maxv <- max(lfc)
halfv <- 0
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.35),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.35))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "snow1", "snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))


pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk)

ctx <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/cetuximab/pdo_cetuxi.tsv', sep='\t', header=T, row.names = 1)
annot_cols <- ctx[, 1, drop=F]
#pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk, annotation_col = annot_cols)
annot_cols <- annot_cols[order(annot_cols$CTG_5000),, drop=F]
annot_cols <- annot_cols[rownames(annot_cols) %in% colnames(lfc),, drop=F]
lfc <- lfc[, match(rownames(annot_cols), colnames(lfc))]

pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk, annotation_col = annot_cols)

####
samples_data$n <- make.unique(paste0(samples_data$sample,samples_data$treat))
te <- t(expr)
mmm <- merge(samples_data, te, by='row.names')
mmm[, colnames(mmm)=="FGFR1" || colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="FGFR2" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="ATOh1" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="OGT" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="STAT3" | colnames(mmm) == "n", drop=F]

# tmm
expr2 <- read.table(gzfile('tmm.tsv.gz'), sep="\t", header=T)
rownames(expr2) <- substr(rownames(expr2), 3, nchar(rownames(expr2)))
te <- t(expr2)
mmm <- merge(samples_data, te, by='row.names')
mmm[, colnames(mmm)=="FGFR1" || colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="FGFR2" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="ATOh1" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="OGT" | colnames(mmm) == "n", drop=F]
mmm[, colnames(mmm)=="STAT3" | colnames(mmm) == "n", drop=F]


deg_expr <- expr2[rownames(expr2) %in% rownames(selection),]
samples_data <- samples_data[order(samples_data$sample, samples_data$treat),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
annot_row <- data.frame(row.names=rownames(deg), expr=log(deg$baseMean), log2FoldChange=deg$log2FoldChange)
annot_row <- annot_row[rownames(annot_row) %in% rownames(deg_expr),]
pheatmap(log(deg_expr+PC), show_rownames=TRUE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)


gene_fc <- function(sample, expr, samples_data) {
  meta <- samples_data[samples_data$sample == sample,]
  expr <- expr[,colnames(expr) %in% rownames(meta)]
  ctx <- meta[meta$treat =="cetuxi",]
  nt <- meta[meta$treat =="NT",]
  ave_ctx <- rowMeans(expr[,colnames(expr) %in% rownames(ctx)])
  ave_nt <- rowMeans(expr[,colnames(expr) %in% rownames(nt)])
  return(ave_ctx/ave_nt)
}

us <- unique(samples_data$sample)

samplefc <- sapply(us, gene_fc, deg_expr, samples_data)
colnames(samplefc) <- us
pheatmap(samplefc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F)

lfc <- log2(samplefc)
minv <- min(lfc)
maxv <- max(lfc)
halfv <- 0
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.35),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.35))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "snow1", "snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))


pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk)

ctx <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/cetuximab/pdo_cetuxi.tsv', sep='\t', header=T, row.names = 1)
annot_cols <- ctx[, 1, drop=F]
#pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk, annotation_col = annot_cols)
annot_cols <- annot_cols[order(annot_cols$CTG_5000),, drop=F]
annot_cols <- annot_cols[rownames(annot_cols) %in% colnames(lfc),, drop=F]
lfc <- lfc[, match(rownames(annot_cols), colnames(lfc))]

pheatmap(lfc, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=F, cluster_rows = F, color = my_palette, breaks = bk, annotation_col = annot_cols)

