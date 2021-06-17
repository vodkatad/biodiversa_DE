library(clusterProfiler)
data <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_ATRX_PDO/ATRX_cutoff0.05-MUT.vs.WT.gseain.tsv', sep="\t", header=TRUE)
geneList <- data$sort
names(geneList) <- as.character(data$geneid)
geneList <- sort(geneList, decreasing = TRUE)
library(msigdbr)


m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, human_gene_symbol)

#minSize=15, maxSize=1000, nperm=10000,
em2 <- GSEA(geneList, TERM2GENE = m_t2g, minGSSize=15, maxGSSize=1000, nPerm=10000, pvalueCutoff=0.1, verbose=FALSE)

library(enrichplot)
gseaplot2(em2, geneSetID = 1, title = em2$Description[1])

# TODO decide n based on sign?
#  p.adjust   qvalues 
#n <- sum(em2$p.adjust < 0.05)
n <- 30
ridgeplot(em2, showCategory=20)
dotplot(em2, showCategory=20)
gseaplot2(em2, geneSetID = 10, title = em2$Description[10])

