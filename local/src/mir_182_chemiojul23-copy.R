library(multiMiR)

p3 <- get_multimir(mirna = 'hsa-miR-182-3p', summary = TRUE)
p3 <- p3@data
#p3 <- p3[c(1:50),]

p5 <- get_multimir(mirna = 'hsa-miR-182-5p', summary = TRUE)
p5 <- p5@data
#p5 <- p5[c(1:272),]

mirna <- unique(c(p3$target_symbol, p5$target_symbol))

deg <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/vsd.tsv.gz"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(deg) <- gsub("H_", "", rownames(deg))
deg$genes <- rownames(deg)
deg <- deg %>% filter(genes %in% mirna)
deg$genes <- NULL

samples <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
samples <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples$genealogy <- rownames(samples)

deg <- as.data.frame(t(deg))
#deg$sample <- substr(rownames(deg), 1 , 7)
deg$genealogy <- rownames(deg)
deg <- merge(deg, samples, by = "genealogy")
deg$batch <- NULL
deg$sample <- NULL
rownames(deg) <- deg$genealogy
deg <- deg[order(deg$type),]
annot <- deg[,c(1, 210)]
annot$genealogy <- NULL
deg$genealogy <- NULL
deg$type <- NULL
deg <- as.data.frame(t(deg))

pheatmap(deg, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annot, show_colnames = FALSE, show_rownames = FALSE)

type <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
type <- read.table(type, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(type) <- gsub("H_", "", rownames(type))
type$genes <- rownames(type)

type$istarget <- ifelse(rownames(type) %in% mirna, 'yes', 'no')

phyper(nrow(type[type$istarget=="yes" & type$padj < 0.05,])-1, nrow(type[type$padj < 0.05,]), nrow(type)-nrow(type[type$padj < 0.05,]), nrow(type[type$istarget=="yes",]), lower.tail=F)

library(ggplot2)
ggplot(data=type, aes(x=log2FoldChange, y=-log10(padj), color=istarget))+geom_point(size=1)+theme_bw()+
  scale_color_manual(values=c('black', 'red'))+geom_hline(yintercept=-log10(0.05))+geom_vline(xintercept=log2(1.5))+
  geom_vline(xintercept=-log2(1.5))




type <- type %>% filter(genes %in% mirna)
type <- type %>% filter(padj < 0.05)

