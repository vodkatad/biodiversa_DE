
load('/scratch/trcanmed/DE_RNASeq/dataset/connector_A_B/methy.Rdata')
merged2 <- merged2[order(merged2$col),]
data <- merged2[, seq(1, (ncol(merged2)-1))]
annot <- data.frame(row.names=rownames(merged2), group=merged2$col)

pheatmap(data, annotation_row = annot, cluster_rows =  T, show_rownames = FALSE, show_colnames = FALSE)

krt <- data[, grepl('KRT', colnames(data))]
pheatmap(krt, annotation_row = annot, cluster_rows =  F, show_rownames = FALSE, show_colnames = FALSE)

load('/scratch/trcanmed/DE_RNASeq/dataset/connector_A_B/dds_B.Rdata')
plotCophenetic()