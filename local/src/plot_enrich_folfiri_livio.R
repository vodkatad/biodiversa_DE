gene_res_f <- "/scratch/trcanmed/DE_RNASeq/dataset/deg_folfiri_new/type_cutoff0.05-PD.vs.PR.gseain_2.tsv"

gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)
type <- "H"

m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici

em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)
#save.image("gsea_results.R")
sign <- as.data.frame(em$Description)
sign$order <- 1:nrow(sign) 

livio <- read.xlsx("/scratch/trcanmed/DE_RNASeq/dataset/deg_folfiri_new/signature_livio.xlsx")
livio <- livio$Signature

sign <- sign %>% filter(`em$Description` %in% livio)

plot_list <- list()

for (i in sign$order) {
  p <- gseaplot2(em, geneSetID = i, title = em$Description[i])
  plot_list[[i]] <- p
  }

pdf("/home/mferri/plots_enrich.pdf")
for (i in 1:17) {
  print(plot_list[[i]])
}
dev.off()

type <- "C2"

m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici

emc2 <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)
signc2 <- as.data.frame(emc2$Description)
signc2$order <- 1:nrow(signc2) 
signc2 <- signc2 %>% filter(`emc2$Description` %in% livio)

plot_list_c2 <- list()

for (i in signc2$order) {
  p <- gseaplot2(emc2, geneSetID = i, title = emc2$Description[i])
  plot_list_c2[[i]] <- p
}

pdf("/home/mferri/plots_enrich_c2.pdf")
for (i in 1:847) {
  print(plot_list_c2[[i]])
}
dev.off()


type <- "C6"

m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici

emc6 <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)
signc6 <- as.data.frame(emc6$Description)
signc6$order <- 1:nrow(signc6) 
signc6 <- signc6 %>% filter(`emc6$Description` %in% livio)

plot_list_c6 <- list()

for (i in signc6$order) {
  p <- gseaplot2(emc6, geneSetID = i, title = emc6$Description[i])
  plot_list_c6[[i]] <- p
}

pdf("/home/mferri/plots_enrich_c6.pdf")
for (i in 1:25) {
  print(plot_list_c6[[i]])
}
dev.off()


