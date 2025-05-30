
library(ggplot2)
CRC0322 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/322_creni_vs_EGF/trattamento_cutoff0.05-Creni_5uM_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/327_creni_vs_EGF/trattamento_cutoff0.05-Creni_5uM_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327$sign <- ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange > 0.58, 'up', ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange < -0.58,'down', 'grigi'))
m <- merge(CRC0327, CRC0322, by='row.names')

ggplot(data=m, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_smooth(method="lm")+geom_point(aes(color=sign))


CRC0322 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/322_combo_vs_EGF/trattamento_cutoff0.05-Combo_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/327_combo_vs_EGF/trattamento_cutoff0.05-Combo_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327$sign <- ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange > 0.58, 'up', ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange < -0.58,'down', 'grigi'))
m <- merge(CRC0327, CRC0322, by='row.names')

ggplot(data=m, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_smooth(method="lm")+geom_point(aes(color=sign))


CRC0322 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/322_CTX_vs_EGF/trattamento_cutoff0.05-Cetux_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327 <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/327_CTX_vs_EGF/trattamento_cutoff0.05-Cetux_72h.vs.EGF0.1.deseq2.tsv', sep="\t", header=T)


CRC0327$sign <- ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange > 0.58, 'up', ifelse(CRC0327$padj < 0.05 & CRC0327$log2FoldChange < -0.58,'down', 'grigi'))
m <- merge(CRC0327, CRC0322, by='row.names')

ggplot(data=m, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_smooth(method="lm")+geom_point(aes(color=sign))
