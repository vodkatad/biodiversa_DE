library(tidyverse)
library(ggplot2)
library(cowplot)

cris_preproc <- snakemake@input[["cris_f"]]
deg_f <- snakemake@input[["samples_f"]]
cris_result_pdf <- snakemake@output[["pdf"]]
log_f <- snakemake@log[['log']]

#cris_final <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_model_cris-right_validated.tsv"
#cris_final <- read.table(cris_final, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#cris_preproc <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv"
cris <- read.table(cris_preproc, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#cris <- as.data.frame(cbind(cris$model, cris$NA..7))
#colnames(cris) <- c("model", "CRIS_PDX")
#cris <- cris %>% distinct()

#deg_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
deg <- read.table(deg_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi <- deg$sample

cris <- cris %>% filter(genealogy %in% casi)

names(cris)[colnames(cris) == "genealogy"] <- "sample"
#deg <- na.omit(deg)

merged <- merge(deg, cris, by="sample")
merged$batch <- NULL
non_responder <- merged %>% filter(type == "non_responder_3Q")
responder <- merged %>% filter(type == "responder_1Q")
non_responder <- non_responder %>% filter(!cris == "HET")
responder <- responder %>% filter(!cris == "HET")


nr <- ggplot(non_responder, aes(x=cris, fill=cris))+geom_bar()+ggtitle("non_responder_3Q")
r <- ggplot(responder, aes(x=cris, fill=cris))+geom_bar()+ggtitle("responder_1Q")

colors <- c('darkorange1', '#d2232de6', '#1d295be6', '#128549e6', '#19ac9be6')
names(colors) <- c('CRIS-A', 'CRIS-B', 'CRIS-C', 'CRIS-D', 'CRIS-E')
legend <- get_legend(nr)+theme(legend.position = "right")
merged <- merged %>% filter(!cris=="HET")
#plot_grid(r + theme(legend.position = "none"), nr, legend)
ggplot(merged, aes(x=cris, fill=cris))+geom_bar()+facet_wrap(~factor(type, levels=c('non_responder_3Q', "responder_1Q")))+scale_fill_manual(values=colors)
ggsave(file=cris_result_pdf)

table_responder <- as.data.frame(table(responder$cris))
table_non_responder <- as.data.frame(table(non_responder$cris))
df <- cbind(table_responder, table_non_responder)
colnames(df) <- c("CRIS", "Freq_res", "CRIS", "Freq_non_res")
df$CRIS <- NULL
rownames(df) <- df$CRIS
df$CRIS <- NULL

sink(log_f, append=TRUE)
crisa <- df %>% filter(rownames(df) == "CRIS-A")
#df <- data.frame(responder=table_responder$Freq, non_responder=table_non_responder$Freq)
crisa <- rbind(crisa, c(31-2, 31-9))
rownames(crisa) <- c("CRISA", "Other")
fisher.test(crisa)
chisq.test(crisa)

crisb <- df %>% filter(rownames(df) == "CRIS-B")
crisb <- rbind(crisb, c(31-7, 31-2))
rownames(crisb) <- c("CRISB", "Other")
fisher.test(crisb)
chisq.test(crisb)

crisc <- df %>% filter(rownames(df) == "CRIS-C")
crisc <- rbind(crisc, c(31-11, 31-10))
rownames(crisc) <- c("CRISC", "Other")
fisher.test(crisc)
chisq.test(crisc)

crisd <- df %>% filter(rownames(df) == "CRIS-D")
crisd <- rbind(crisd, c(31-8, 31-2))
rownames(crisd) <- c("CRISD", "Other")
fisher.test(crisd)
chisq.test(crisc)

crise <- df %>% filter(rownames(df) == "CRIS-E")
crise <- rbind(crise, c(31-3, 31-8))
rownames(crise) <- c("CRISE", "Other")
fisher.test(crise)
chisq.test(crisd)
sink()



