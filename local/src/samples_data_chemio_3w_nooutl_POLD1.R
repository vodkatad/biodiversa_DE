library(tidyverse)
library(DESeq2)
library(ggplot2)

dds <- snakemake@input[["dds_originale"]]
samples <- snakemake@input[["samples_3w"]]
meta <- snakemake@output[["sample"]]

#dds <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/dds.Rdata"
load(dds)
what <- "type"

data<-plotCounts(dds, "H_POLD1", intgroup=what, returnData=T)
data <- data[order(data$count),]
e2 <- ggplot(data, aes_string(x = what, y = "count"))
e3 <- e2 + geom_jitter(aes_string(shape = what, color = what),   position = position_jitter(0.2),size = 3) + stat_summary( aes_string(color = what), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.4, color="darkgreen")+theme_bw()+scale_y_continuous(trans='log10')+labs(color = "Xeno", x="Xeno", shape="Xeno", y="Log10(nreads)")+ggtitle("POLD1")
e4 <- e3 + geom_text(aes(label = rownames(data)), vjust = -0.5, hjust = 0.5, size = 3)
#ggsave(paste0("POLD1.eps"), width = 150, height = 100, units = "mm")

data <- data %>% filter(count < 288.6516)

pold1_outliers <- rownames(data)

#samples <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
samples <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples$genealogy <- rownames(samples)
samples <- samples %>% filter(!genealogy %in% pold1_outliers)
samples$genealogy <- NULL

write.table(samples, file=meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)