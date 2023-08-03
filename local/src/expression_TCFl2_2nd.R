fpkmtot <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/fpkm.tsv.gz", quote = "", sep = "\t",
                      header = TRUE, stringsAsFactors = FALSE)
geni <- c("ATOH1", "KRT8", "GAPDH", "PCDHGB2", "MREG", "MTPAP", "ATP5E")

fpkmtot <- fpkmtot %>% filter(rownames(fpkmtot) %in% geni)

# fpkm80 <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/CRC0080/CRC0080_fpkm.tsv.gz", quote = "", sep = "\t",
#                      header = TRUE, stringsAsFactors = FALSE)
# 
# fpkm80 <- fpkm80 %>% filter(rownames(fpkm80) %in% geni)
# 
# c <- cor.test(fpkm80$CRC0080_NE_R1, fpkmtot$CRC0080_NE_R1)

fpkmtot <- as.data.frame(t(fpkmtot))
fpkmtot$case <- rownames(fpkmtot)
fpkmtot <- fpkmtot[order(fpkmtot$case),]
fpkmtot$geno <- substr(fpkmtot$case, 9, 10)
fpkmtotne <- fpkmtot %>% filter(geno == "NE")
fpkmtotn2 <- fpkmtot %>% filter(geno=="N2")

fpkmtotfin <- rbind(fpkmtotne, fpkmtotn2)
fpkmtotfin$case <- NULL
fpkmtotfin$geno <- NULL

fpkmtot$case <- NULL
fpkmtot$geno <- NULL
pheatmap(log2(fpkmtotfin+1), cluster_rows = FALSE, cluster_cols = FALSE,fontsize_row = 5)
