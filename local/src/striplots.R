#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)

input <- args[[1]]
ddsf <- args[[2]]
outputplotdir <- args[[3]]
what <- args[[4]]


load(ddsf)

striplots <- function(gene, dds) {
    if (!any(rownames(dds) == gene)) {
        return()
     }
      data<-plotCounts(dds, gene, intgroup=what, returnData=T)
    data <- data[order(data$count),]
    e2 <- ggplot(data, aes_string(x = what, y = "count"))
      e3 <- e2 + geom_jitter(  aes_string(shape = what, color = what),   position = position_jitter(0.2),size = 3) + stat_summary( aes_string(color = what), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.4, color="darkgreen")+theme_bw()+scale_y_continuous(trans='log10')+labs(color = "Xeno", x="Xeno", shape="Xeno", y="Log10(nreads)")
      ggsave(paste0(gene, ".png"))
}

#run on philae, (rnaseq) /home/grassi/RNAseq_biodiversa/dataset/DESeq

list <- read.table(input, header=F, sep="\t")
list$V1 <- paste0("H_", list$V1)
setwd(outputplotdir)
garbage <- lapply(list$V1, striplots, dds)
