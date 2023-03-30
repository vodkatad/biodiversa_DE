load('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/magnum_deg/H_GSEA.Rdata')
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

n <- nrow(em[em@result$p.adjust < 0.01,])

ems <- em[em@result$p.adjust< 0.01,]
table(ems$enrichmentScore < 0)

#dotplot(em, showCategory=n, color="pvalue", x="GeneRatio", font.size=10)ridgeplot(em, showCategory = 20, fill="pvalue", core_enrichment = TRUE)
ridgeplot(em, showCategory = n, fill="pvalue", core_enrichment = TRUE)

ws <-c('HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_HEDGEHOG_SIGNALING', 'HALLMARK_ANGIOGENESIS')
wi <- c()
for (w in ws) {
  i <- which(em$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}
gseaplot2(em, geneSetID=wi, subplots = 1:2)
##
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)

load('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/magnum_deg/C2_GSEA.Rdata')

n <- nrow(em[em@result$p.adjust < 0.05,])

#dotplot(em, showCategory=n, color="pvalue", x="GeneRatio", font.size=10)ridgeplot(em, showCategory = 20, fill="pvalue", core_enrichment = TRUE)
ridgeplot(em, showCategory = 20, fill="pvalue", core_enrichment = TRUE)

react <- em[grepl('REACTOME', em@result$ID),, asis=T]
n <- nrow(react[react@result$p.adjust < 0.05,])
ridgeplot(react, showCategory = 20, fill="pvalue", core_enrichment = TRUE)


sreact <- react[react@result$p.adjust < 0.05,, asis=T]

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/magnum/all_matrix_reactome_scores.tsv', sep="\t", header=T)

show <- d[rownames(d)%in% sreact$Description, ]


annot <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/magnum_deg/samples_data', sep="\t", header=T, row.names=1)
annot <- annot[, c('response'), drop=F]
annot <- annot[order(annot$response), , drop=F]
show <- show[, match(rownames(annot), colnames(show))]
pheatmap(as.matrix(show), annotation_col = annot, show_colnames = F, cluster_cols = F, show_rownames = F)


sreact <- em[em@result$p.adjust < 0.05,, asis=T]
sreact <- sreact[grepl('REACTOME', sreact@result$ID),, asis=T]

sreact <- sreact[order(-sreact$NES),]
sreactdown <- tail(sreact, n=20)
sreactup <- sreact[sreact$NES > 0,]
sreactup <- sreactup[order(-sreactup$NES),]
sreact <- rbind(sreactdown, sreactup)
show <- d[rownames(d)%in% sreact$Description, ]

show <- show[, match(rownames(annot), colnames(show))]
pheatmap(as.matrix(show), annotation_col = annot, show_colnames = F, cluster_cols = F, show_rownames = T)

ws <-c('REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE', 'REACTOME_KERATINIZATION')
wi <- c()
for (w in ws) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}
gseaplot2(em, geneSetID=wi, subplots = 1:2)asd