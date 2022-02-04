library(limma)
library(tidyverse)
setwd("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584/")

load("paperoga.Rdata")
efilterGenes <- rowSums(new_data > minc) < minsamples
edata <- new_data[!efilterGenes,]
e_new_data <- edata[,match(rownames(metadata), colnames(edata))]

model <- model.matrix(as.formula("~batch+type+0"), data=metadata)
png("voom.png")
voomnorm <- voom(e_new_data, model, plot=T)
dev.off()
fit <- lmFit(voomnorm, model)
contr <- makeContrasts(typeLMO_BASALE - typeLMX_BASALE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
write.table(top.table, "limma.tsv", sep="\t", col.names=T, quote=F)

lmolmx <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584/type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE.deseq2.tsv"
ox <- read.table(lmolmx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

signlimma <- rownames(top.table[top.table$adj.P.Val<0.05,])
signox <- rownames(ox[ox$padj<0.05,])

inters <- length(intersect(singlimma, signox))

merged <- merge(ox, top.table, by = "row.names")
 
plot(merged$log2FoldChange, merged$logFC)

merged_genes <- as.data.frame(merged[merged$log2FoldChange < -4,]$Row.names)
colnames(merged_genes) <- c("genes")
merged_genes <- merged_genes %>% mutate(genes = gsub("H_", "", genes))

write.table(merged_genes, file = "genes_limma.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

