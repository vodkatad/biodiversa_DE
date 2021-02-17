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

#THR <- 0.5849625
THR <- 1
PC <- 0.0001
deg_lf <- deg[deg$log2FoldChange > THR,]

deg_expr <- expr[rownames(expr) %in% rownames(deg_lf),]

samples_data <- samples_data[order(samples_data$sample, samples_data$treat),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
pheatmap(log(deg_expr+PC), show_rownames=TRUE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F)
annot_row <- data.frame(row.names=rownames(deg), expr=log(deg$baseMean), log2FoldChange=deg$log2FoldChange)

annot_row <- annot_row[rownames(annot_row) %in% rownames(deg_expr),]
pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)

annot_row <- annot_row[annot_row$expr>4.791111,]
deg_expr <- deg_expr[rownames(deg_expr) %in% rownames(annot_row),]
pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)

#### merge
deg_up <- deg[deg$log2FoldChange > 0,]
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
write.table(res, gzfile('big_merge.tsv.gz'), sep="\t", quote=F, row.names=T)

test <- res[(res$log2FC_PDO_R > -0.01 & res$log2FC_PDO_R < 0.01 & res$padj_basal_NT > 0.05) | is.na(res$log2FC_PDO_R) | is.na(res$padj_basal_NT),]

deg_expr <- expr[rownames(expr) %in% rownames(test),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)


test <- res[res$padj_notpaired < 0.05 | is.na(res$padj_notpaired),]

deg_expr <- expr[rownames(expr) %in% rownames(test),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)



test <- res[(res$log2FC_PDO_R > -0.01 & res$log2FC_PDO_R < 0.01 & res$padj_basal_NT > 0.05 & res$padj_notpaired < 0.05 & !is.na(res$padj_notpaired)) | is.na(res$log2FC_PDO_R) | is.na(res$padj_basal_NT),]

deg_expr <- expr[rownames(expr) %in% rownames(test),]
deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)



