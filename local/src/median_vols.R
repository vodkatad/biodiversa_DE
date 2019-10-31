#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]
outputplot <- args[[3]]


indf <- read.table(input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
indf <- indf[!is.na(indf$vol), ]

terziles <- quantile(indf$vol, c(1/3, 2/3))
terziles
#median <- median(indf$vol)
#median
indf$chemio <- ifelse(indf$vol <= terziles[1], "resp", ifelse(indf$vol >= terziles[2], "nonresp", "grey"))
indf$chemio <- factor(indf$chemio, levels=c("nonresp","grey","resp"))
ggplot(indf, aes(y=vol,x=reorder(id, -vol),fill=chemio))+geom_col()+ylab("Delta %volume")+xlab("Case")+theme_bw()+theme(axis.text.x = element_blank())+scale_fill_manual(values=c("red","grey","blue"))+ggtitle(paste0("Thrs =",paste0(terziles,collapse=", ")))
ggsave(outputplot)

table(indf$chemio)
indf$basale <- NULL
indf$type <- NULL
indf$vol <- NULL
table(indf$chemio, indf$batch)
table(indf$chemio, indf$coverage)


write.table(indf, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
