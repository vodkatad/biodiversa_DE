### AB vs all per krt significative

krt <- "/scratch/trcanmed/DE_RNASeq/dataset/connector_Aa_vsall/validated_keratin.tsv"
krt <- read.table(krt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
krt <- unique(krt$SYMBOL)

abba <- "/scratch/trcanmed/DE_RNASeq/dataset/connector_Ab_allPD_revision/Ab_Ba/Cluster_cutoff0.05-Ab.vs.Ba.goinsplit_down_Ab_Ba.tsv"
abba <- read.table(abba, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
abbau <- "/scratch/trcanmed/DE_RNASeq/dataset/connector_Ab_allPD_revision/Ab_Ba/Cluster_cutoff0.05-Ab.vs.Ba.goinsplit_up_Ab_Ba.tsv"
abbau <- read.table(abbau, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

ba <- "/scratch/trcanmed/DE_RNASeq/dataset/connector_Ab_allPD_revision/Ab_Ba/Cluster_cutoff0.05-Ab.vs.Ba.deseq2_Ab_Ba.tsv"
ba <- read.table(ba, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(ba)
rownames(ba) <- NULL
ba <- cbind(genes,ba)
ba <- ba %>% mutate(genes = gsub("H_", "", genes))
ba <- ba %>% filter(log2FoldChange < 0)
ba <- ba %>% filter(padj < 0.05)

ba_krt <- ba
ba_krt <- ba_krt %>% filter(genes %in% krt)
ba_krt <- ba_krt$genes

bb <- "/scratch/trcanmed/DE_RNASeq/dataset/connector_Ab_allPD_revision/Ab_Bb/Cluster_cutoff0.05-Ab.vs.Bb.deseq2_Ab_Bb.tsv"
bb <- read.table(bb, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(bb)
rownames(bb) <- NULL
bb <- cbind(genes,bb)
bb <- bb %>% mutate(genes = gsub("H_", "", genes))
bb <- bb %>% filter(log2FoldChange < 0)
bb <- bb %>% filter(padj < 0.05)

bb_krt <- bb
bb_krt <- bb_krt %>% filter(genes %in% krt)
bb_krt <- bb_krt$genes

diffbabb <- setdiff(ba_krt, bb_krt)
diffbbba <- setdiff(bb_krt, ba_krt)
intersection <- intersect(ba_krt, bb_krt)

krt_sign <- data.frame(matrix(ncol = 3, nrow = 34))
colnames(krt_sign) <- c("genes", "sign_Ab.vs.Ba", "sign_Ab.vs.Bb")
krt_sign$genes <- krt

for (i in krt_sign$genes) {
  if (i %in% ba_krt) {
    krt_sign[krt_sign$genes == i, "sign_Ab.vs.Ba"] <- TRUE
  } else if (!i %in% ba_krt) {
    krt_sign[krt_sign$genes == i, "sign_Ab.vs.Ba"] <- FALSE
  } 
}

for (i in krt_sign$genes) {
  if (i %in% bb_krt) {
    krt_sign[krt_sign$genes == i, "sign_Ab.vs.Bb"] <- TRUE
  } else if (!i %in% bb_krt) {
    krt_sign[krt_sign$genes == i, "sign_Ab.vs.Bb"] <- FALSE
  } 
}

krt_sign$both <- NA

for (i in krt_sign$genes) {
  if (krt_sign[krt_sign$genes == i, "sign_Ab.vs.Ba"] == TRUE & krt_sign[krt_sign$genes == i, "sign_Ab.vs.Bb"] == TRUE) {
    krt_sign[krt_sign$genes == i, "both"] <- TRUE
  } else {
    krt_sign[krt_sign$genes == i, "both"] <- FALSE
  }
}

write.table(krt_sign, file = "/home/mferri/significant_krt_Ab_revision.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
