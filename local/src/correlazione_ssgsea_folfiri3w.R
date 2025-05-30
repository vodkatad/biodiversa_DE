## ssgsea correlazione folfiri

library(ggsignif)

old <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_H_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
old <- read.table(old, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
old <- old %>% filter(p.adjust < 0.05)

gsea <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/h-scores_tot_fpkm.tsv"
gsea <- read.table(gsea, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea <- as.data.frame(t(gsea))
gsea$sample <- rownames(gsea)

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

merged <- merge(gsea, res, by="sample")
rownames(merged) <- merged$sample
merged$sample <- NULL

results <- data.frame(
  name = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:50) {
  cor_test <- cor.test(merged[[i]], merged$`3wks`, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  results <- rbind(results, data.frame(
    name = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

diz <- merged
diz$`3wks` <- NULL
diz <- as.data.frame(t(diz))
diz$gsea <- rownames(diz)
diz$name <- 1:nrow(diz)
diz <- diz[,c(92,93)]

results <- merge(results, diz, by="name")
results$name <- NULL 
results <- results[,c(3,1,2)]

forplot <- merged
forplot$genealogy <- rownames(forplot)

samplesdata <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data", quote = "",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samplesdata$genealogy <- rownames(samplesdata)
samplesdata <- samplesdata[,c("genealogy", "type")]

forplot <- merge(forplot, samplesdata, by="genealogy", all.x=TRUE)
rownames(forplot) <- forplot$genealogy
forplot$genealogy <- NULL
forplot$type[is.na(forplot$type)] <- "no_deg"
forplot$type <- as.factor(forplot$type)

forwaterfall <- forplot
forwaterfall$labels <- rownames(forwaterfall)

ggplot(forwaterfall, aes(x = reorder(labels, -HALLMARK_IL6_JAK_STAT3_SIGNALING), y = HALLMARK_IL6_JAK_STAT3_SIGNALING, fill=type)) +
  geom_point(stat = "identity", aes(color=type)) + scale_color_manual(values = c("grey", "red", "blue"))

ggplot(forwaterfall, aes(x = reorder(labels, -`3wks`), y = HALLMARK_INTERFERON_ALPHA_RESPONSE, fill=type)) +
  geom_point(stat = "identity", aes(color=type)) + scale_color_manual(values = c("grey", "red", "blue"))

ggplot(forwaterfall, aes(x = reorder(labels, -`3wks`), y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, fill=type)) +
  geom_point(stat = "identity", aes(color=type)) + scale_color_manual(values = c("grey", "red", "blue"))


# HALLMARK_IL6_JAK_STAT_SIGNALING
# HALLMARK_INTERFERON_ALPHA_RESPONSE
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_EPITHELIAL_MESENCHIMAL_TRANSITION
# HALLMARK_FATTY_ACID_METABOLISM
# HALLMARK_KRAS_SIGNALING_UP

for (i in rownames(forplot)) {
  if (forplot[i,"type"]=="responder_1Q") {
    forplot[i,"type_2"] <- "Sensitive_1Q"
  } else if (forplot[i, "type"] == "non_responder_3Q") {
    forplot[i, "type_2"] <- "Resistant_3Q"
  } else {
    forplot[i,"type_2"] <- "Intermediate_Response_2Q"
  }
}

forplot$type <- forplot$type_2
forplot$type_2 <- NULL

setwd("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/")
pdf("HALLMARK_IL6_JAK_STAT3_SIGNALING.pdf")
p <- ggplot(data=forplot, aes(x=type, y=HALLMARK_IL6_JAK_STAT3_SIGNALING))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p,"HALLMARK_IL6_JAK_STAT3_SIGNALING.rds" )
dev.off()
pdf("HALLMARK_INTERFERON_ALPHA_RESPONSE.pdf")
p <- ggplot(data=forplot, aes(x=type, y=HALLMARK_INTERFERON_ALPHA_RESPONSE))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "HALLMARK_INTERFERON_ALPHA_RESPONSE.rds")
dev.off()
pdf("HALLMARK_INTERFERON_GAMMA_RESPONSE.pdf")
p <- ggplot(data=forplot, aes(x=type, y=HALLMARK_INTERFERON_GAMMA_RESPONSE))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p,"HALLMARK_INTERFERON_GAMMA_RESPONSE.rds" )
dev.off()
pdf("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf")
p <- ggplot(data=forplot, aes(x=type, y=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.rds")
dev.off()
pdf("HALLMARK_FATTY_ACID_METABOLISM.pdf")
p <- ggplot(data=forplot, aes(x=type, y=HALLMARK_FATTY_ACID_METABOLISM))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "HALLMARK_FATTY_ACID_METABOLISM.rds")
dev.off()
#ggplot(forplot, aes(x=HALLMARK_IL6_JAK_STAT3_SIGNALING, y=`3wks`))+geom_point(aes(color=type))+geom_smooth(method = "lm")
#ggplot(forplot, aes(x=HALLMARK_APICAL_JUNCTION, y=`3wks`))+geom_point(aes(color=type))+geom_smooth(method = "lm")

## c2

gsea <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c2-scores_tot_fpkm.tsv"
gsea <- read.table(gsea, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea <- as.data.frame(t(gsea))
gsea$sample <- rownames(gsea)

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

merged <- merge(gsea, res, by="sample")
rownames(merged) <- merged$sample
merged$sample <- NULL

results <- data.frame(
  name = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

merged <- as.data.frame(t(merged))
merged$sign <- rownames(merged)
signature <- c("REACTOME_FATTY_ACID_METABOLISM", "LIN_APC_TARGETS", "KIM_MYC_AMPLIFICATION_TARGETS_DN", "SANSOM_APC_TARGETS_DN", "3wks")
merged <- merged %>% filter(sign %in% signature)
merged$sign <- NULL
merged <- as.data.frame(t(merged))

for (i in 1:4) {
  cor_test <- cor.test(merged[[i]], merged$`3wks`, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  results <- rbind(results, data.frame(
    name = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

diz <- merged
diz$`3wks` <- NULL
diz <- as.data.frame(t(diz))
diz$gsea <- rownames(diz)
diz$name <- 1:nrow(diz)
diz <- diz[,c(92,93)]

results <- merge(results, diz, by="name")
results$name <- NULL 
results <- results[,c(3,1,2)]

forplot <- merged
forplot$genealogy <- rownames(forplot)

samplesdata <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data", quote = "",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samplesdata$genealogy <- rownames(samplesdata)
samplesdata <- samplesdata[,c("genealogy", "type")]

forplot <- merge(forplot, samplesdata, by="genealogy", all.x=TRUE)
rownames(forplot) <- forplot$genealogy
forplot$genealogy <- NULL
forplot$type[is.na(forplot$type)] <- "no_deg"
forplot$type <- as.factor(forplot$type)

# "REACTOME_FATTY_ACID_METABOLISM", "LIN_APC_TARGETS", "KIM_MYC_AMPLIFICATION_TARGETS_DN", "SANSOM_APC_TARGETS_DN"
for (i in rownames(forplot)) {
  if (forplot[i,"type"]=="responder_1Q") {
    forplot[i,"type_2"] <- "Sensitive_1Q"
  } else if (forplot[i, "type"] == "non_responder_3Q") {
    forplot[i, "type_2"] <- "Resistant_3Q"
  } else {
    forplot[i,"type_2"] <- "Intermediate_Response_2Q"
  }
}

forplot$type <- forplot$type_2
forplot$type_2 <- NULL

pdf("REACTOME_FATTY_ACID_METABOLISM.pdf")
p <- ggplot(data=forplot, aes(x=type, y=REACTOME_FATTY_ACID_METABOLISM))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "REACTOME_FATTY_ACID_METABOLISM.rds")
dev.off()
pdf("LIN_APC_TARGETS.pdf")
p <- ggplot(data=forplot, aes(x=type, y=LIN_APC_TARGETS))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "LIN_APC_TARGETS.rds")
dev.off()
pdf("KIM_MYC_AMPLIFICATION_TARGETS_DN.pdf")
p <- ggplot(data=forplot, aes(x=type, y=KIM_MYC_AMPLIFICATION_TARGETS_DN))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "KIM_MYC_AMPLIFICATION_TARGETS_DN.rds")
dev.off()
pdf("SANSOM_APC_TARGETS_DN.pdf")
p <- ggplot(data=forplot, aes(x=type, y=SANSOM_APC_TARGETS_DN))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "SANSOM_APC_TARGETS_DN.rds")
dev.off()
## c6

gsea <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c6-scores_tot_fpkm.tsv"
gsea <- read.table(gsea, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea <- as.data.frame(t(gsea))
gsea$sample <- rownames(gsea)

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

merged <- merge(gsea, res, by="sample")
rownames(merged) <- merged$sample
merged$sample <- NULL

results <- data.frame(
  name = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

merged <- as.data.frame(t(merged))
merged$sign <- rownames(merged)
signature <- c("BCAT.100_UP.V1_UP", "LEF1_UP.V1_UP", "LEF1_UP.V1_DN", "3wks")
merged <- merged %>% filter(sign %in% signature)
merged$sign <- NULL
merged <- as.data.frame(t(merged))

for (i in 1:3) {
  cor_test <- cor.test(merged[[i]], merged$`3wks`, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  results <- rbind(results, data.frame(
    name = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

diz <- merged
diz$`3wks` <- NULL
diz <- as.data.frame(t(diz))
diz$gsea <- rownames(diz)
diz$name <- 1:nrow(diz)
diz <- diz[,c(92,93)]

results <- merge(results, diz, by="name")
results$name <- NULL 
results <- results[,c(3,1,2)]

forplot <- merged
forplot$genealogy <- rownames(forplot)

samplesdata <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data", quote = "",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samplesdata$genealogy <- rownames(samplesdata)
samplesdata <- samplesdata[,c("genealogy", "type")]

forplot <- merge(forplot, samplesdata, by="genealogy", all.x=TRUE)
rownames(forplot) <- forplot$genealogy
forplot$genealogy <- NULL
forplot$type[is.na(forplot$type)] <- "no_deg"
forplot$type <- as.factor(forplot$type)

for (i in rownames(forplot)) {
  if (forplot[i,"type"]=="responder_1Q") {
    forplot[i,"type_2"] <- "Sensitive_1Q"
  } else if (forplot[i, "type"] == "non_responder_3Q") {
    forplot[i, "type_2"] <- "Resistant_3Q"
  } else {
    forplot[i,"type_2"] <- "Intermediate_Response_2Q"
  }
}
forplot$type <- forplot$type_2
forplot$type_2 <- NULL
# "BCAT.100_UP.V1_UP", "LEF1_UP.V1_UP", "LEF1_UP.V1_DN"
pdf("BCAT.100_UP.V1_UP.pdf")
p <- ggplot(data=forplot, aes(x=type, y=BCAT.100_UP.V1_UP))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "BCAT.100_UP.V1_UP.rds")
dev.off()
pdf("LEF1_UP.V1_UP.pdf")
p <- ggplot(data=forplot, aes(x=type, y=LEF1_UP.V1_UP))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "LEF1_UP.V1_UP.rds")
dev.off()
pdf("LEF1_UP.V1_DN.pdf")
p <- ggplot(data=forplot, aes(x=type, y=LEF1_UP.V1_DN))+geom_boxplot(aes(col = type))+geom_signif(comparisons = list(c("Resistant_3Q", "Sensitive_1Q")))+scale_color_manual(values = c("grey", rgb(165, 0, 25, maxColorValue = 255), rgb(30, 85, 130, maxColorValue = 255)))
saveRDS(p, "LEF1_UP.V1_DN.rds")
dev.off()
