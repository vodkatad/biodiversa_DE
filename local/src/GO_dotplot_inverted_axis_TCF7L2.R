library(dplyr)
library(ggplot2)

## MUT_N2.vs.NE
# s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_samples_data",
#                 quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# length(unique(s$model))
# 
# deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
#                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# deg$genes <- rownames(deg)
# deg2 <- deg %>% filter(padj < 0.05)
# deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
# 
# up <- deg %>% filter(padj<0.05)
# up <- up %>% filter(log2FoldChange > 0.5849625)
# 
# down <- deg %>% filter(padj < 0.05)
# down <- down %>% filter(log2FoldChange < -0.5849625)


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(y="GeneRatio", x="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  xlab(NULL) +
  ylim(0.01, 0.2)+
  scale_size(range=c(3, 8))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('~/TCF7L2_mut_N2_vs_NE_GO_down.pdf')

## WT_N2.vs.NE
# s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_samples_data",
#                 quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# length(unique(s$model))
# 
# deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
#                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# deg$genes <- rownames(deg)
# deg2 <- deg %>% filter(padj < 0.05)
# deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
# up <- deg %>% filter(padj<0.05)
# up <- up %>% filter(log2FoldChange > 0.5849625)
# 
# down <- deg %>% filter(padj < 0.05)
# down <- down %>% filter(log2FoldChange < -0.5849625)

go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))


ggplot(df, aes_string(y="GeneRatio", x="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  xlab(NULL) +
  ylim(0.01, 0.2)+
  scale_size(range=c(3, 8))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('~/TCF7L2_WT_N2_vs_NE_GO_down.pdf')
