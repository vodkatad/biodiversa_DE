library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/onf_pdo/pp_cutoff0.05-plastic.vs.rigid.gseain_2.tsv"

gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)
type <- "C2"

# prendere sta cazzo di DNA_METHYLATION a mano TODO

m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, gene_symbol) ### altrimenti chiede gli id numerici

em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)
#save.image("gsea_results.R")
sign <- as.data.frame(em$Description)
sign$order <- 1:nrow(sign) 

ss <- em@result
w <- c('REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 'REACTOME_METHYLATION', 'REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP')
wi <- c()
for (s in w) {
  wi <- c(wi, grep(s, rownames(ss)))
}

plot_list <- list()

j <- 1
for (i in wi) {
  p <- gseaplot2(em, geneSetID = i, title = em$Description[i])
  print(em$Description[i])
  print(em$NES[i])
  print(em$qvalue[i])
  plot_list[[j]] <- p
  j <- j+1
}

###########
type <- "C3"

# prendere sta cazzo di DNA_METHYLATION a mano TODO

m_t2g <- msigdbr(species = "Homo sapiens", category = type) %>% 
  dplyr::select(gs_name, gene_symbol) ### altrimenti chiede gli id numerici

em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)
#save.image("gsea_results.R")
sign <- as.data.frame(em$Description)
sign$order <- 1:nrow(sign) 


dd <- as.data.frame(em@result)

p <- gseaplot2(em, geneSetID = 9, title = em$Description[9])

#### let's get methy signature from v2026
make_gene_set <- function(file) {
  sign <- read.table(file, sep="\t", header=F)
  sign$V1
}

make_gene_sets <- function(files, names) {
  res <- lapply(files, make_gene_set)
  names(res) <- names
  geneIdsToGeneSetCollection(res)
}

make_gene_set_fgsea <- function(file, name) {
  sign <- read.table(file, sep="\t", header=F)
  data.frame(term=name, name=sign$V1)
}

#new_signs <- make_gene_sets(c('~/epigenetic_regulation.tsv', '~/dna_methylation.tsv'), c('REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 'REACTOME_DNA_METHYLATION'))
new_signs <- mapply(make_gene_set_fgsea, c('~/epigenetic_regulation.tsv', '~/dna_methylation.tsv'), c('REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 'REACTOME_DNA_METHYLATION'), SIMPLIFY=FALSE, USE.NAMES=FALSE)
new_signs_df <- do.call(rbind, new_signs)

em <- GSEA(geneList, TERM2GENE = new_signs_df, pvalueCutoff = 1, minGSSize=1, maxGSSize=1000)
plot_list <- list()
for (i in 1:length(new_signs)) {
  p <- gseaplot2(em, geneSetID = i, title = em$Description[i])
  print(em$Description[i])
  print(em$NES[i])
  print(em$qvalue[i])
  plot_list[[i]] <- p
}


