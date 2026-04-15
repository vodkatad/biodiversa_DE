ne <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/geno_mut_cutoff0.05-NE_MUT.vs.NE_WT.deseq2.tsv"
ne <- read.table(ne, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ne$genes <- rownames(ne)
ne <- ne[,c("log2FoldChange", "genes", "padj")]
names(ne)[names(ne)=="log2FoldChange"] <- "log2FoldChange_NE"
names(ne)[names(ne)=="padj"] <- "padj_NE"

## lmo like tcf

bb <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv"
bb <- read.table(bb, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bb$genes <- gsub("H_", "", rownames(bb))
bb <- bb[,c("log2FoldChange", "genes", "padj")]
names(bb)[names(bb)=="log2FoldChange"] <- "log2FoldChange_biob_LMO"
names(bb)[names(bb)=="padj"] <- "padj_biob_LMO"

merged <- merge(ne, bb, by="genes")
merged$significance <- NA  

for (i in 1:nrow(merged)) {
  if (merged$padj_NE[i] < 0.05 & merged$padj_biob_LMO[i] < 0.05) {
    merged$significance[i] <- "both"
  } else if (merged$padj_NE[i] < 0.05) {
    merged$significance[i] <- "NE"
  } else if (merged$padj_biob_LMO[i] < 0.05) {
    merged$significance[i] <- "biob_LMO"
  } else {
    merged$significance[i] <- "no_sign"
  }
}

cor.test(merged$log2FoldChange_NE, merged$log2FoldChange_biob_LMO)

ggplot(merged, aes(x = log2FoldChange_NE, 
                   y = log2FoldChange_biob_LMO, 
                   color = significance)) +
  geom_point() +
  geom_smooth(method = 'lm', size = 1, se = FALSE, color = "black") +
  scale_color_manual(values = c("both"     = "green",
                                "NE"       = "blue",
                                "biob_LMO" = "red",
                                "no_sign"  = "lightgrey"))


## lmx like TCF

bb <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF_LMX/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv"
bb <- read.table(bb, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bb$genes <- gsub("H_", "", rownames(bb))
bb <- bb[,c("log2FoldChange", "genes", "padj")]
names(bb)[names(bb)=="log2FoldChange"] <- "log2FoldChange_biob_LMX"
names(bb)[names(bb)=="padj"] <- "padj_biob_LMX"

merged <- merge(ne, bb, by="genes")
merged$significance <- NA  

for (i in 1:nrow(merged)) {
  if (merged$padj_NE[i] < 0.05 & merged$padj_biob_LMX[i] < 0.05) {
    merged$significance[i] <- "both"
  } else if (merged$padj_NE[i] < 0.05) {
    merged$significance[i] <- "NE"
  } else if (merged$padj_biob_LMX[i] < 0.05) {
    merged$significance[i] <- "biob_LMX"
  } else {
    merged$significance[i] <- "no_sign"
  }
}

cor.test(merged$log2FoldChange_NE, merged$log2FoldChange_biob_LMX)

ggplot(merged, aes(x = log2FoldChange_NE, 
                   y = log2FoldChange_biob_LMX, 
                   color = significance)) +
  geom_point() +
  geom_smooth(method = 'lm', size = 1, se = FALSE, color = "black") +
  scale_color_manual(values = c("both"     = "green",
                                "NE"       = "blue",
                                "biob_LMX" = "red",
                                "no_sign"  = "lightgrey"))

## extended

bb <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_MUT.vs.WT_biobanca/geno_cutoff0.05-MUT.vs.WT.deseq2.tsv"
bb <- read.table(bb, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bb$genes <- gsub("H_", "", rownames(bb))
bb <- bb[,c("log2FoldChange", "genes", "padj")]
names(bb)[names(bb)=="log2FoldChange"] <- "log2FoldChange_biob_extended_LMX"
names(bb)[names(bb)=="padj"] <- "padj_biob_extended_LMX"

merged <- merge(ne, bb, by="genes")
merged$significance <- NA  

for (i in 1:nrow(merged)) {
  if (merged$padj_NE[i] < 0.05 & merged$padj_biob_extended_LMX[i] < 0.05) {
    merged$significance[i] <- "both"
  } else if (merged$padj_NE[i] < 0.05) {
    merged$significance[i] <- "NE"
  } else if (merged$padj_biob_extended_LMX[i] < 0.05) {
    merged$significance[i] <- "biob_ext_LMX"
  } else {
    merged$significance[i] <- "no_sign"
  }
}

cor.test(merged$log2FoldChange_NE, merged$log2FoldChange_biob_extended_LMX)

ggplot(merged, aes(x = log2FoldChange_NE, 
                   y = log2FoldChange_biob_extended_LMX, 
                   color = significance)) +
  geom_point() +
  geom_smooth(method = 'lm', size = 1, se = FALSE, color = "black") +
  scale_color_manual(values = c("both"     = "green",
                                "NE"       = "blue",
                                "biob_ext_LMX" = "red",
                                "no_sign"  = "lightgrey"))

