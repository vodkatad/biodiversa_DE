old <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
old <- read.table(old, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
old <- old %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5849625)
#diff <- gsub("H_", "", rownames(old))
up <- old %>% filter(log2FoldChange > 0.5849625)
down <- old %>% filter(log2FoldChange < - 0.5849625)
diff_up <- gsub("H_", "", rownames(up))
diff_down <- gsub("H_", "", rownames(down))

casi_f <- "/scratch/trcanmed/DE_RNASeq/local/share/data/chemio_def_jul23/CHEMIO_WATERFALL_PLOT_Eugy_Luglio2023.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi$ff30 <- NULL
casi <- casi[order(casi$X3WKS, decreasing = TRUE),]

metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMX_BASALE", type))
meda_f$CASE <- substr(meda_f$sample_id_R, 1,7)
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
merged <- merge(casi, meda_f, by = "CASE")
merged <- merged[,c("sample_id_R", "X3WKS", "CASE")]
res <- merged
colnames(res) <- c("sample", "3wks", "CASE")
res <- res %>% filter(!CASE == "CRC0578")
res$CASE <- NULL

casi <- res$sample

vsd <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/vsd.tsv.gz"

vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$genes <- gsub("H_", "", rownames(vsd))
vsd_up <- vsd %>% filter(genes %in% diff_up)
rownames(vsd_up) <- vsd_up$genes
vsd_up$genes <- NULL
vsd_up <- as.data.frame(t(vsd_up))
vsd_up$sample <- rownames(vsd_up)
vsd_up <- merge(vsd_up, res, by="sample")
rownames(vsd_up) <- vsd_up$sample
vsd_up$sample <- NULL


vsd_down <- vsd %>% filter(genes %in% diff_down)
rownames(vsd_down) <- vsd_down$genes
vsd_down$genes <- NULL
vsd_down <- as.data.frame(t(vsd_down))
vsd_down$sample <- rownames(vsd_down)
vsd_down <- merge(vsd_down, res, by="sample")
rownames(vsd_down) <- vsd_down$sample
vsd_down$sample <- NULL

results_up <- data.frame(
  name = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

results_down <- data.frame(
  name = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)


for (i in colnames(vsd_up)[1:22]) {
  cor_test <- cor.test(vsd_up[[i]], vsd_up$`3wks`, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  results_up <- rbind(results_up, data.frame(
    name = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

for (i in colnames(vsd_down)[1:40]) {
  cor_test <- cor.test(vsd_down[[i]], vsd_down$`3wks`, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  results_down <- rbind(results_down, data.frame(
    name = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

for (i in rownames(results_up)) {
  if (results_up[i, "p_value"] < 0.05) {
    results_up[i,"color"] <- "yes"
  } else {
    results_up[i, "color"] <- "no"
  }
}

for (i in rownames(results_down)) {
  if (results_down[i, "p_value"] < 0.05) {
    results_down[i, "color"] <- "yes"
  } else {
    results_down[i, "color"] <- "no"
  }
}

ggplot(results_up, aes(x=estimate))+geom_histogram()
ggplot(results_down, aes(x=estimate))+geom_histogram()

res <- rbind(results_up, results_down)
res$name <- factor(res$name, levels = res$name[order(-res$estimate)])

ggplot(res, aes(x=name, y=estimate, color=color))+ geom_jitter(height=0, shape=18)+scale_color_manual(values =c("black", "red"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(results_up, aes(x=name, y=estimate, color=color))+ geom_jitter(height=0, shape=18)+scale_color_manual(values =c("black", "red"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(results_down, aes(x=name, y=estimate, color=color))+geom_jitter(height = 0, shape=18)+scale_color_manual(values = c("black", "red"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


results_up_sing <- results_up %>% filter(p_value < 0.05)
results_down_sign <- results_down %>% filter(p_value < 0.05)
