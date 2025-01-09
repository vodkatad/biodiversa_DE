

mary <- read.table('/mnt/trcanmed/snaketree/stash/NEW_DESeq2_Response_3wkPD_v_PR_at_Placebo.csv', sep=",", header=T, row.names=1)
us <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2_geni_scelti.tsv', sep="\t", header=T, row.names = 1)

m <- merge(mary, us, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL

ggplot(data=m, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_point()+geom_smooth(method="lm")+theme_bw()



mary6 <- read.table('/mnt/trcanmed/snaketree/stash/NEW_DESeq2_Response_6wkPD_v_PR_at_Placebo.csv', sep=",", header=T, row.names=1)

m <- merge(mary6, us, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL

ggplot(data=m, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_point()+geom_smooth(method="lm")+theme_bw()
cor.test(m$log2FoldChange.x, m$log2FoldChange.y)
