## correlazione nuovi grigi - old magnifici 29

new <- "/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
new <- read.table(new, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
new$genes <- gsub("H_", "", rownames(new))
names(new)[colnames(new)=="log2FoldChange"] <- "log2FoldChange_m29_grigi"
names(new)[colnames(new)=="padj"] <- "padj_m29_grigi"
names(new)[colnames(new)=="pvalue"] <- "pvalue_m29_grigi"

#old <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
#old <- read.table(old, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#old$genes <- gsub("H_", "", rownames(old))

old <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
old <- read.table(old, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
old$genes <- gsub("H_", "", rownames(old))
names(old)[colnames(old)=="log2FoldChange"] <- "log2FoldChange_m29"
names(old)[colnames(old)=="padj"] <- "padj_m29"
names(old)[colnames(old)=="pvalue"] <- "pvalue_m29"

merged <- merge(new, old, by="genes")

ci <- cor.test(merged$log2FoldChange_m29_grigi, merged$log2FoldChange_m29)
ci[["estimate"]]
ci[["p.value"]]


sign_m29_grigi <- merged %>% filter(padj_m29_grigi < 0.05)

ggplot(sign_m29_grigi, aes(x=log2FoldChange_m29_grigi, y=log2FoldChange_m29, col = ifelse(padj_m29 < 0.05,'yes','no')))+
   geom_point(size=1)+guides(color=guide_legend(title="Significative in m29_grigi
  significative in m29"))#+geom_text(label=merged$genes, hjust=0, vjust=0)

chegenisono <- sign_m29_grigi
chegenisono <- chegenisono %>% filter(log2FoldChange_m29_grigi > 0)
chegenisono <- chegenisono %>% filter(log2FoldChange_m29 < 0)

sign_m29 <- merged %>% filter(padj_m29 < 0.05)
 
ggplot(sign_m29, aes(x=log2FoldChange_m29_grigi, y=log2FoldChange_m29, col = ifelse(padj_m29_grigi < 0.05,'yes','no')))+
  geom_point(size=1)+guides(color=guide_legend(title="Significativi in m29 
  significative in m29_grigi"))#+geom_text(label=merged$genes, hjust=0, vjust=0)

cm29_grigi <- cor.test(sign_m29$log2FoldChange_m29_grigi, sign_m29$log2FoldChange_m29)
cm29_grigi[["estimate"]]
cm29_grigi[["p.value"]]

rownames(sign_m29) <- sign_m29$genes

for (i in rownames(sign_m29)) {
  if (sign_m29[i, "log2FoldChange_m29"] * sign_m29[i, "log2FoldChange_m29_grigi"]>0) {
    sign_m29[i, "concordanza"] <- "concordi"
  } else {
    sign_m29[i, "concordanza"] <- "discordi"
  }
}

for (i in rownames(sign_m29)) {
  if (sign_m29[i, "padj_m29_grigi"]< 0.05) {
    sign_m29[i, "significativi_grigi"] <- "sì" 
  } else {
    sign_m29[i,  "significativi_grigi"] <- "no"
  }
}

sign_m29 <- merged %>% filter(pvalue_m29 < 0.05)

ggplot(sign_m29, aes(x=log2FoldChange_m29_grigi, y=log2FoldChange_m29, col = ifelse(padj_m29_grigi < 0.05,'yes','no')))+
  geom_point(size=1)+guides(color=guide_legend(title="Significativi in m29 
  significative in m29_grigi"))#+geom_text(label=merged$genes, hjust=0, vjust=0)

cm29_grigi <- cor.test(sign_m29$log2FoldChange_m29_grigi, sign_m29$log2FoldChange_m29)
cm29_grigi[["estimate"]]
cm29_grigi[["p.value"]]

rownames(sign_m29) <- sign_m29$genes

for (i in rownames(sign_m29)) {
  if (sign_m29[i, "log2FoldChange_m29"] * sign_m29[i, "log2FoldChange_m29_grigi"]>0) {
    sign_m29[i, "concordanza"] <- "concordi"
  } else {
    sign_m29[i, "concordanza"] <- "discordi"
  }
}

for (i in rownames(sign_m29)) {
  if (sign_m29[i, "padj_m29_grigi"]< 0.05) {
    sign_m29[i, "significativi_grigi"] <- "sì" 
  } else {
    sign_m29[i,  "significativi_grigi"] <- "no"
  }
}

# remove <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/removefromDEG.tsv")
# 
# h <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv", quote = "",
#                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# c6 <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv", quote = "",
#                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# c2 <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv", quote = "",
#                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# h <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv", quote = "",
#                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# c6 <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv", quote = "",
#                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

w3 <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
w3 <- read.table(w3, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
w3 <- w3 %>% filter(padj < 0.05)
w3 <- w3 %>% filter(abs(log2FoldChange) > 0.5849625)
w3$genes <- gsub("H_", "", rownames(w3))

names(w3)[colnames(w3)=="log2FoldChange"] <- "log2FoldChange_w3"
names(new)[colnames(new)=="log2FoldChange"] <- "log2FoldChange_m29_subgroup"

names(w3)[colnames(w3)=="padj"] <- "padj_3w"
names(w3)[colnames(w3)=="pvalue"] <- "pvalue_3w"

names(new)[colnames(new)=="padj"] <- "padj_m29_subgroup"
names(new)[colnames(new)=="pvalue"] <- "pvalue_m29_subgroup"

## per venn up down
w3_up <- w3 %>% filter(log2FoldChange_w3 > 0.5849625)
w3_down <- w3 %>% filter(log2FoldChange_w3 < -0.5849625)

new_up <- new %>% filter(genes %in% w3$genes)
new_down <- new %>% filter(genes %in% w3$genes)
new_up <- new_up %>% filter(log2FoldChange_m29_subgroup > 0.5849625)
new_down <- new_down %>% filter(log2FoldChange_m29_subgroup < -0.5849625)

venn_up <- intersect(w3_up$genes, new_up$genes)
venn_down <- intersect(w3_down$genes, new_down$genes)

petros <- merge(w3, new, by="genes")
ci <- cor.test(petros$log2FoldChange_w3, petros$log2FoldChange_m29_subgroup)
ci[["estimate"]]
ci[["p.value"]]
rownames(petros) <- petros$genes

for (i in rownames(petros)) {
  if (petros[i, "log2FoldChange_w3"]* petros[i,"log2FoldChange_m29_grigi"] > 0) {
    petros[i, "direction"] <- "concorde"
  } else {
    petros[i, "direction"] <- "discorde"
  }
}

length(setdiff(w3$genes, rownames(petros)))

petros <- petros[,c("log2FoldChange_w3", "pvalue_3w", "padj_3w", "log2FoldChange_m29_grigi", "pvalue_m29_grigi", "padj_m29_grigi", "direction")]

cambia <- petros %>% filter(direction == "discorde")

petros$segno <- NULL
cambia$segno <- NULL


write.xlsx(petros, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/correlazione_petros_m29_subgroup.xlsx", rowNames=TRUE)
write.xlsx(cambia, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/correlazione_petros_m29_subgroup_cambianosegno.xlsx", rowNames=TRUE)
