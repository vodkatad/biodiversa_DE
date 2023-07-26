d <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm_H.tsv.gz'), sep="\t", header=TRUE, stringsAsFactors = FALSE)

means <- apply(d, 1, mean)
medianMean <- median(means)

overMedMea <- means[means > medianMean]

dexpr <- d[names(overMedMea),]

library(dplyr)

data <- as.data.frame(overMedMea)
data$quantile <- ntile(data$overMedMea, n=4)

w <- c('ATOH1','DLL1','GFI1','DEFA5','DEFA6','KRT8','GAPDH','PCDHGB2','MREG','MTPAP','ATP5E','APCDD1','AXIN1')

data[rownames(data)%in% w,]