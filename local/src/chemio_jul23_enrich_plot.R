library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(showtext)
library(RColorBrewer)

load("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/gsea_results_h.R")

h <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_H_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv",
                quote = "",sep = "\t", header = TRUE)

ws <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_MYC_TARGETS_V1","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN")

h <- h %>% filter(ID %in% ws)
hpos <- h %>% filter(enrichmentScore > 0)
hpos <- hpos %>% filter(pvalue < 0.05)
hneg <- h %>% filter(enrichmentScore < 0)
hneg <- hneg %>% filter(pvalue < 0.05)

wspos <- hpos$ID
wsneg <- hneg$ID

wi <- c()

for (w in wspos) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

p <- gseaplot2(em, geneSetID=wi, color = brewer.pal(4, "Set2"), subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/hallmark_positive_chemiojul23_scritte.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = brewer.pal(4, "Set2"), subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/hallmark_positive_chemiojul23_muto.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)

wi <- c()
for (w in wsneg) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

p <- gseaplot2(em, geneSetID=wi, color = c("#A6D854", "#FFD92F"), subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/hallmark_negative_chemiojul23_scritte.jpeg", width=170, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("#A6D854", "#FFD92F"), subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/hallmark_negative_chemiojul23_muto.jpeg", width=170, height=89, units="mm", device="jpeg", dpi=300)


load("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/gsea_results_c2.R")

c2 <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C2_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv",
                quote = "",sep = "\t", header = TRUE)


# > ws <-c('REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE', 'REACTOME_KERATINIZATION')
# > wi <- c()
# > for (w in ws) {
#   +   i <- which(em@result$Description==w)
#   +   wi <- c(wi, i)
#   +   #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
#     + }
# > gseaplot2(em, geneSetID=wi, subplots = 1:2)

ws <- c("WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
        "REACTOME_FATTY_ACID_METABOLISM","KEGG_FATTY_ACID_METABOLISM",
        "WP_TCA_CYCLE_AKA_KREBS_OR_CITRIC_ACID_CYCLE","REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
        "REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION","KIM_MYC_AMPLIFICATION_TARGETS_DN",
        "LIN_APC_TARGETS","SANSOM_APC_TARGETS_DN")

# c2 <- c2 %>% filter(ID %in% ws)
# c2pos <- c2 %>% filter(enrichmentScore > 0)
# c2pos <- c2pos %>% filter(pvalue < 0.05)
# c2neg <- c2 %>% filter(enrichmentScore < 0)
# c2neg <- c2neg %>% filter(pvalue < 0.05)

c2 <- c2 %>% filter(ID %in% ws)
c2 <- c2 %>% filter(pvalue < 0.05)
#wspos <- c2pos$ID
#wsneg <- c2neg$ID
c2 <- c2$ID

wi <- c()

for (w in c2) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

p <- gseaplot2(em, geneSetID=wi, subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c2_chemiojul23_scritte.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c2_chemiojul23_muto.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)

# wi <- c()
# for (w in wsneg) {
#   i <- which(em@result$Description==w)
#   wi <- c(wi, i)
#   #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
# }
# 
# p <- gseaplot2(em, geneSetID=wi, subplots = 1:2)
# ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c2_negative_chemiojul23_scritte.jpeg", width=170, height=89, units="mm", device="jpeg", dpi=300)
# p <- mygseaplot(em, geneSetID=wi, subplots = 1:2)
# ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c2_negative_chemiojul23_muto.jpeg", width=170, height=89, units="mm", device="jpeg", dpi=300)


load("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/gsea_results_c6.R")

c6 <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C6_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv",
                 quote = "",sep = "\t", header = TRUE)
# > ws <-c('REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE', 'REACTOME_KERATINIZATION')
# > wi <- c()
# > for (w in ws) {
#   +   i <- which(em@result$Description==w)
#   +   wi <- c(wi, i)
#   +   #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
#     + }
# > gseaplot2(em, geneSetID=wi, subplots = 1:2)
ws <- c("LEF1_UP.V1_UP","LEF1_UP.V1_DN","BCAT.100_UP.V1_UP")
c6 <- c6 %>% filter(ID %in% ws)
c6 <- c6 %>% filter(pvalue < 0.05)
#wspos <- c2pos$ID
#wsneg <- c2neg$ID
c2 <- c2$ID


wi <- c()

for (w in ws) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

p <- gseaplot2(em, geneSetID=wi, subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c6_chemiojul23_scritte.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, subplots = 1:2)
ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/c6_chemiojul23_muto.jpeg", width=200, height=89, units="mm", device="jpeg", dpi=300)
