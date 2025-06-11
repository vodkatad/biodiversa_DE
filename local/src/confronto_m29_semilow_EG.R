m29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
m29 <- read.table(m29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
m29$genes <- gsub("H_", "", rownames(m29))

rad51 <- "/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
rad51 <- read.table(rad51, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rad51$genes <- gsub("H_", "", rownames(rad51))

merged <- merge(m29, rad51, by="genes")

cor.test(merged$log2FoldChange.x, merged$log2FoldChange.y)

ggplot(data=merged, aes(x=log2FoldChange.x, y=log2FoldChange.y))+geom_point()+geom_smooth()+theme_bw(base_size=20)

h_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv"
h_m29 <- read.table(h_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

h_r_f <- "/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv"
h_r <- read.table(h_r_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# h_w3 <- h_w3 %>% filter(p.adjust < 0.05)
# 
# h_m29 <- h_m29 %>% filter(p.adjust < 0.05)
ws <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_MYC_TARGETS_V1","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN")

h_r <- h_r %>% filter(ID %in% ws)
h_m29 <- h_m29 %>% filter(ID %in% ws)
names(h_r)[names(h_r) == 'enrichmentScore'] <- "enrichmentScore_r"
names(h_m29)[names(h_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"

merged_h <- merge(h_r, h_m29, by="ID")
merged_h <- merged_h[c("ID", "enrichmentScore_r", "enrichmentScore_m29")]

sink(log_f, append=TRUE)
"Correlation enrichment Score gsea H"
cor.test(merged_h$enrichmentScore_r, merged_h$enrichmentScore_m29)
sink()

# Pearson's product-moment correlation
# 
# data:  merged_h$enrichmentScore_w3 and merged_h$enrichmentScore_m29
# t = 2.4594, df = 6, p-value = 0.04916
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.00772091 0.94258852
# sample estimates:
#       cor 
# 0.7085389 

write.table(merged_h, file=cor_gsea_h, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

c2_r_f <- "/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv"
c2_r <- read.table(c2_r_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)  

c2_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv"
c2_m29 <- read.table(c2_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# c2_w3 <- c2_w3 %>% filter(p.adjust < 0.05)
# 
# c2_m29 <- c2_m29 %>% filter(p.adjust < 0.05)

ws <- c("REACTOME_FATTY_ACID_METABOLISM","KIM_MYC_AMPLIFICATION_TARGETS_DN",
        "LIN_APC_TARGETS","SANSOM_APC_TARGETS_DN")

c2_r <- c2_r %>% filter(ID %in% ws)
c2_m29 <- c2_m29 %>% filter(ID %in% ws)
names(c2_r)[names(c2_r) == 'enrichmentScore'] <- "enrichmentScore_r"
names(c2_m29)[names(c2_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"

merged_c2 <- merge(c2_r, c2_m29, by="ID")
merged_c2 <- merged_c2[c("ID", "enrichmentScore_r", "enrichmentScore_m29")]

write.table(merged_c2, file=cor_gsea_c2, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
sink(log_f, append=TRUE)
"Correlation enrichment Score gsea c2"
cor.test(merged_c2$enrichmentScore_r, merged_c2$enrichmentScore_m29)
sink()
# Pearson's product-moment correlation
# 
# data:  merged_c2$enrichmentScore_w3 and merged_c2$enrichmentScore_m29
# t = 2.4642, df = 7, p-value = 0.0432
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.03184686 0.92636879
# sample estimates:
#       cor 
# 0.6815536


c6_r_f <- "/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv"
c6_r <- read.table(c6_r_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)  

c6_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv"
c6_m29 <- read.table(c6_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#c6_w3 <- c6_w3 %>% filter(p.adjust < 0.05)

#c6_m29 <- c6_m29 %>% filter(p.adjust < 0.05)
ws <- c("LEF1_UP.V1_UP","LEF1_UP.V1_DN","BCAT.100_UP.V1_UP")

c6_r <- c6_r %>% filter(ID %in% ws)
c6_m29 <- c6_m29 %>% filter(ID %in% ws)
names(c6_r)[names(c6_r) == 'enrichmentScore'] <- "enrichmentScore_r"
names(c6_m29)[names(c6_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"

merged_c6 <- merge(c6_r, c6_m29, by="ID")
merged_c6 <- merged_c6[c("ID", "enrichmentScore_r", "enrichmentScore_m29")]

write.table(merged_c6, file=cor_gsea_c6, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
sink(log_f, append=TRUE)
"Correlation enrichment Score gsea C6"
cor.test(merged_c6$enrichmentScore_r, merged_c6$enrichmentScore_m29)
sink()

fra <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/IHC_m29_merged.tsv"
fra <- read.table(fra, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
fra$model <- rownames(fra)

sm29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/samples_data"
sm29 <- read.table(sm29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sm29$batch <- NULL
sm29 <- sm29 %>% filter(!duplicated(model))

m29_fra <- merge(fra, sm29, by="model")

m29_fra_res <- m29_fra %>% filter(type=="resistant")
media_rad51_resistenti_m29 <- mean(m29_fra_res$RAD51_NT)

m29_fra_sens <- m29_fra %>% filter(type=="sensitive")
media_rad51_sensitive_m29 <- mean(m29_fra_sens$RAD51_NT)

m29_fra_nohighrad <- m29_fra %>% filter(!model %in% c("CRC0029", "CRC0151", "CRC0479", "CRC0204"))

m29_fra_res_nohighrad <- m29_fra_nohighrad %>% filter(type=="resistant")
media_rad51_resistenti_semilow <- mean(m29_fra_res_nohighrad$RAD51_NT)

m29_fra_sens_nohighrad <- m29_fra_nohighrad %>% filter(type=="sensitive")
media_rad51_sensitive_semilow <- mean(m29_fra_sens_nohighrad$RAD51_NT)

t.test(m29_fra_res_nohighrad$RAD51_NT, m29_fra_sens_nohighrad$RAD51_NT, alternative = "two.sided", var.equal = FALSE)

ggplot(m29_fra_nohighrad, aes(x=type, y=RAD51_NT))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw(base_size=20)

### rimozione CRC0077
m29_fra_nohighrad <- m29_fra %>% filter(!model %in% c("CRC0077","CRC0029", "CRC0151", "CRC0479", "CRC0204"))

m29_fra_res_nohighrad <- m29_fra_nohighrad %>% filter(type=="resistant")
media_rad51_resistenti_semilow <- mean(m29_fra_res_nohighrad$RAD51_NT)

m29_fra_sens_nohighrad <- m29_fra_nohighrad %>% filter(type=="sensitive")
media_rad51_sensitive_semilow <- mean(m29_fra_sens_nohighrad$RAD51_NT)

t.test(m29_fra_res_nohighrad$RAD51_NT, m29_fra_sens_nohighrad$RAD51_NT, alternative = "two.sided", var.equal = FALSE)

ggplot(m29_fra_nohighrad, aes(x=type, y=RAD51_NT))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw(base_size=20)

## venn bruti a mano, hypg
# egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ wc -l *up* *down*
#   33 type_cutoff0.05-resistant.vs.sensitive.goinsplit_up.tsv
# 11 type_cutoff0.05-resistant.vs.sensitive.goinsplit_down.tsv
# 44 total
# egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ wc -l ../magnifici29/*up*
#   wc: ../magnifici29/barplot_go_type_cutoff0.05-resistant.vs.sensitive_up: Is a directory
# 0 ../magnifici29/barplot_go_type_cutoff0.05-resistant.vs.sensitive_up
# 761 ../magnifici29/GO_results_type_cutoff0.05-resistant.vs.sensitive_up.tsv
# 38 ../magnifici29/type_cutoff0.05-resistant.vs.sensitive.goinsplit_up.tsv
# 799 total
# egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ wc -l ../magnifici29/*down*
#   wc: ../magnifici29/barplot_go_type_cutoff0.05-resistant.vs.sensitive_down: Is a directory
# 0 ../magnifici29/barplot_go_type_cutoff0.05-resistant.vs.sensitive_down
# 790 ../magnifici29/GO_results_type_cutoff0.05-resistant.vs.sensitive_down.tsv
# 13 ../magnifici29/type_cutoff0.05-resistant.vs.sensitive.goinsplit_down.tsv
# 803 total
# egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ filter_1col 1 type_cutoff0.05-resistant.vs.sensitive.goinsplit_up.tsv < ../magnifici29/type_cutoff0.05-resistant.vs.sensitive.goinsplit_up.tsv
# | wc -l
# 24


#egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ cat ../magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv | filter_1col 1 <(cut -f 1 type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv | sed 1d) | wc -l
#17892

# up, urn is semilow
phyper(23, 33, 17892-33, 38, lower.tail=F)

# down, urn is semilow
phyper(4, 11, 17892-11, 13, lower.tail=F)

#egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/rad51_res.vs.sens$ filter_1col 1 type_cutoff0.05-resistant.vs.sensitive.goinsplit_down.tsv < ../magnifici29/type_cutoff0.05-resistant.vs.sensitive.goinsplit_down.tsv | wc -l
#5


