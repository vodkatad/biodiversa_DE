#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)

input <- args[[1]]
ddsf <- args[[2]]
outputplotdir <- args[[3]]

load(ddsf)

striplots <- function(gene, dds) {
    if (!any(rownames(dds) == gene)) {
        return()
     }
      data<-plotCounts(dds, gene, intgroup="chemio", returnData=T)
  data$chemio <- factor(data$chemio, levels=c("nonresp", "grey", "resp"))
    data <- data[order(data$count),]
    e2 <- ggplot(data, aes(x = chemio, y = count))
      e3 <- e2 + geom_jitter(  aes(shape = chemio, color = chemio),   position = position_jitter(0.2),size = 3) + stat_summary( aes(color = chemio), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.4, color="darkgreen")+     scale_color_manual(values =  c("#FC4E07", "#778899","#00AFBB"))+theme_bw()+scale_y_continuous(trans='log10')+labs(color = "Xeno", x="Xeno", shape="Xeno", y="Log10(nreads)")
      ggsave(paste0(gene, ".png"))
}

#run on philae, (rnaseq) /home/grassi/RNAseq_biodiversa/dataset/DESeq

list <- read.table(input, header=F, sep="\t")
list$V1 <- paste0("H_", list$V1)
setwd(outputplotdir)
garbage <- lapply(list$V1, striplots, dds)
