#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)

input <- args[[1]]
ddsf <- args[[2]]
outputplotdir <- args[[3]]
what <- args[[4]]
samples <- args[[5]]

load(ddsf)
if (samples == 'all') {
   samples <- unique(metadata$sample)
} else {
   samples <- unlist(strsplit(samples, ','))
}
# FIXME change XENO to something dictated by conf.sk 

striplots <- function(gene, dds) {
    if (!any(rownames(dds) == gene)) {
        return()
     }
    data <-plotCounts(dds, gene, intgroup=what, returnData=T, transform=FALSE) # temp transform to false for investigations with Alberto
    gene <- gsub('H_','', gene)
    data$sample <- substr(rownames(data),0,7)
    data <- data[data$sample %in% samples,]
    data <- data[order(data$sample, data$treat, rownames(data)),] # order is preserved like the exp design here? yup!
    data$group <- make.unique(paste0(data$treat, '_', data$sample))
    data$group <- sapply(data$group, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
    e2 <- ggplot(data, aes_string(x = 'treat', y = "count"))
    e2 + geom_point(  aes_string(shape = what, color = 'sample'),  size = 3) +theme_bw()+labs(color = "Xeno", x="Xeno", shape="Xeno", y="nreads")+geom_line(aes_string(x="treat",y="count", color="sample", group="group"))+ggtitle(gene)
    ggsave(paste0(gene, ".png"))
    foldchange <- function(group, data)  {
        d <- data[data$group==group,, drop=F]
        return(d[d$treat=="cetuxi", 'count']/ d[d$treat=="NT", 'count'])
    }

    res <- sapply(unique(data$group), foldchange, data)
    write.table(res, paste0(gene, '.tsv'), sep="\t", quote=F)

}

list <- read.table(input, header=F, sep="\t")
list$V1 <- paste0("H_", list$V1)
setwd(outputplotdir)
garbage <- lapply(list$V1, striplots, dds)
