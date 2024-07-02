setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')

tfeb <- read.table('LMX_BASALE-TFEB_fpkm_ave.tsv', sep="\t", header=TRUE)
nuak <- read.table('LMX_BASALE-NUAK2_fpkm_ave.tsv', sep="\t", header=TRUE)
ulk <- read.table('LMX_BASALE-ULK1_fpkm_ave.tsv', sep="\t", header=TRUE)

m <- merge(tfeb, nuak, by="model")
m1 <- merge(m, ulk, by="model")

colnames(m1) <- c('model', 'tfeb', 'nuak', 'ulk')

cor.test(m1$tfeb, m1$nuak, method='spearman')
cor.test(m1$tfeb, m1$ulk, method='spearman')
cor.test(m1$nuak, m1$ulk, method='spearman')


plot(m1$tfeb, m1$nuak)
plot(m1$tfeb, m1$ulk)

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/fpkm.tsv.gz', sep="\t", header=TRUE)
cor.test(unlist(d['H_NUAK2',]), unlist(d['H_ULK1',]))
plot(unlist(d['H_NUAK2',]), unlist(d['H_ULK1',]))
cor.test(unlist(d['H_NUAK2',]), unlist(d['H_TFEB',]))
cor.test(unlist(d['H_ULK1',]), unlist(d['H_TFEB',]))
plot(unlist(d['H_ULK1',]), unlist(d['H_TFEB',]))


meta <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/samples_data', sep="\t", header=TRUE, row.names=1)
td <- t(d)
m <- merge(td, meta, by="row.names")
library(ggplot2)
ggplot(data=m, aes(x=H_ULK1, y=H_TFEB, color=sample, shape=treat))+geom_point()+theme_bw(base_size=15)

lmfit <- lm(data=m, formula=as.formula("H_ULK1~H_TFEB+sample"))
