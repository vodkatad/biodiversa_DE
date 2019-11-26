#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]]
outputplot <- args[[3]]


indf <- read.table(input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
indf <- indf[!is.na(indf$vol), ]

# TODO get quantiles as params and generalize at the rule level
thrs <- quantile(indf$vol, c(1/8, 7/8))
thrs
#median <- median(indf$vol)
#median
indf$chemio <- ifelse(indf$vol <= thrs[1], "resp", ifelse(indf$vol >= thrs[2], "nonresp", "grey"))
indf$chemio <- factor(indf$chemio, levels=c("nonresp","grey","resp"))
ggplot(indf, aes(y=vol,x=reorder(id, -vol),fill=chemio))+geom_col()+ylab("Delta %volume")+xlab("Case")+theme_bw()+theme(axis.text.x = element_blank())+scale_fill_manual(values=c("red","grey","blue"))+ggtitle(paste0("Thrs =",paste0(thrs,collapse=", ")))
ggsave(outputplot)

table(indf$chemio)
indf$basale <- NULL
indf$type <- NULL
indf$vol <- NULL
table(indf$chemio, indf$batch)
table(indf$chemio, indf$coverage)


write.table(indf, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
