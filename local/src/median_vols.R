#!/usr/bin/env Rscript
library(ggplot2)
options(error=traceback)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]
outputplot <- args[[3]]


indf <- read.table(input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
indf <- indf[!is.na(indf$vol), ]
median <- median(indf$vol)
median
indf$chemio <- ifelse(indf$vol >= median, "nonresp", "resp")
ggplot(indf, aes(y=vol,x=reorder(id, -vol),fill=chemio))+geom_col()+ylab("Delta %volume")+xlab("Case")+theme_bw()+theme(axis.text.x = element_blank())+scale_fill_manual(values=c("red","blue"))+ggtitle(paste0("Median =",median))
ggsave(outputplot)

table(indf$chemio)
indf$basale <- NULL
indf$type <- NULL
indf$vol <- NULL
table(indf$chemio, indf$batch)
table(indf$chemio, indf$coverage)

ggplot(indf, aes(y=vol,x=reorder(id, -vol),fill=RNAseq))+geom_col()+ylab("Delta %volume")+xlab("Case")+theme_bw()+theme(axis.text.x = element_blank())+scale_fill_manual(values=c("gold","black"))

write.table(indf, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
