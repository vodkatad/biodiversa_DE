#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
gene <- args[[2]]
samples <- args[[3]]
outputplothisto <- args[[4]]
outputplotdist <- args[[5]]
outputstats <- args[[6]]
thr <- args[[7]]

if (thr=="yes") {
  thrp <- TRUE 
} else {
  thrp <- FALSE
}

data <- read.table(gzfile(input), sep="\t", header=TRUE, comment.char='')
# select only wanted samples
if (samples != 'all') {
  data <- data[, grepl(samples, colnames(data))]
}
# get only expressed genes
if (thrp) {
  medians <- apply(data, 1, median)
  thr <- median(medians)
  data <- data[medians > thr,]
}

means <- apply(data, 1, mean)
medians <- apply(data, 1, median)
sdss <- apply(data, 1, sd)
allstats <- data.frame(mean=means, medians=medians, sd=sdss)
sums <- do.call(cbind, lapply(allstats, summary))

if (!gene %in% rownames(data)) {
  stop(paste0('cannot find gene ', gene))
}
geneexpr <- as.numeric(data[gene == rownames(data),])
gmean <- mean(geneexpr)
gmedian <- median(geneexpr)
gsd <- sd(geneexpr)

nsamples <- ncol(data)
cat(nsamples)
genedf <- data.frame(expr=geneexpr)
ggplot(genedf, aes(expr)) + geom_histogram()+theme_bw()
ggsave(outputplothisto)

sums <- rbind(sums, c(gmean, gmedian,gsd))
rownames(sums)[nrow(sums)] <- gene
write.table(sums, file=outputstats, col.names=TRUE, row.names=TRUE, sep="\t", quote=F)

long <- melt(allstats)

plot <- function(all, gene, title) {
  minx <- min(all$value)
  maxx <- quantile(all$value, c(99/100))
  if (gene > maxx) {
    maxx <- gene
  }
  return(ggplot(all, aes(x=value)) + geom_density() +xlim(c(minx, maxx)) +geom_vline(xintercept = gene, color = "red")+ggtitle(title))
}
save.image('pippo.RData')
pdf(outputplotdist,onefile = TRUE)
print(plot(long[long$variable=="mean",], gmean, 'mean'))
print(plot(long[long$variable=="medians",], gmedian,'median'))
print(plot(long[long$variable=="sd",], gsd, 'sd'))
dev.off()
