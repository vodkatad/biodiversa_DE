#pca1 = prcomp(object, scale. = TRUE)

## script to change genes list and corrisponding pca
library(DESeq2)

stroma <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/stromal_genes.tsv"
sg <- read.table(stroma, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
vsd_f <- "/scratch/trcanmed//DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# genes <- rownames(vsd)
# rownames(vsd) <- NULL
# vsd <- cbind(genes,vsd)
# vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
# vsd <- vsd %>% remove_rownames %>% column_to_rownames(var="genes")

load("dds.Rdata")
## assay works with DESeq2 and dds.Rdata already loaded!!!
tvsd <- t(assay(vsd)) # è il nostro object
# questa si fa solo la pca

pca_complicata <- function (matrix, ntop = 500){
  cv <- colVars(matrix)
  select <- order(cv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(cv)))]
  matrix <- matrix[, select]
  pca <- prcomp(matrix)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  return(list(pca = pca$x, pca_loading = pca$rotation, percentVar = percentVar))
}

list_pca <- pca_complicata(tvsd)



pca_cool <- function (object, pca, percentVar, intgroup = "type", pc1=1, pc2 =2, returnData = FALSE, plotfile) {
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca[, pc1], PC2 = pca[, pc2], group = group, 
                  intgroup.df, name = colnames(object))
  colnames(d)[1] <- paste0("PC", pc1)
  colnames(d)[2] <- paste0("PC", pc2)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc1:pc2]
    return(d)
  }
  # alpha and size
  p1 <- ggplot(data = d, aes_string(x = colnames(d)[1], y = colnames(d)[2], color = "group")) + 
    geom_point(size = 2, alpha = 0.3) + xlab(paste0(colnames(d)[1],": ", round(percentVar[pc1] * 100), "% variance")) + 
    ylab(paste0(colnames(d)[2],": ", round(percentVar[pc2] * 100), "% variance"))+ coord_fixed() + theme_bw()
  
  #ggsave(filename = plotfile, plot = p1) 
  print(p1)
}
vsd <- vst(dds, blind = FALSE)

pca_cool(pca = list_pca$pca, percentVar = list_pca$percentVar, object = vsd)

#pca_genes <- as.data.frame(list_pca$pca_loading)
#pca_genes <- pca_genes[order(-pca_genes$PC1),]
#delete_genes <- gsub("H_", "", head(rownames(pca_genes), 100))
#write.table(delete_genes, file = "100_geni_pca.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

#stromal_genes <- sg$V1
#tvsd2 <- tvsd[,!colnames(tvsd) %in% delete_genes]
tvsd2 <- tvsd[,!colnames(tvsd) %in% sg$V1]

pca_no_genes <- pca_complicata(tvsd2)
pca_cool(object = vsd, pca = pca_no_genes$pca, percentVar = pca_no_genes$percentVar)

intersect <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/inteserct_stromal_universe.tsv"
uni_claudio <- read.table(intersect, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)





#delete_genes is used both for GO and GSEA and heatmap at te bottom of the code


#geneList <- delete_genes$`head(pca_genes$PC1, 100)`

#for gsea we need also the loading values (as a frequence) to make gsea works
# geneList <- delete_genes
# geneList <- sort(geneList, decreasing = TRUE)

#names(geneList) <- as.character(delete_genes$`gsub("H_", "", head(rownames(pca_genes), 100))`)
#geneList <- sort(geneList, decreasing = TRUE)

geneList <- delete_genes
geneUni <- gsub("H_", "", rownames(list_pca$pca_loading))

egocc <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "CC",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

egomf <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "MF",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

egobp <- enrichGO(gene          = geneList,
                  universe      = geneUni,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",  
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

#dir.create(out_dir)
#setwd(out_dir)

setwd("/home/mferri/")

barplot(egocc, showCategory = 20)
ggsave("CC_100_pca.pdf")
barplot(egomf, showCategory = 20)
ggsave("MF_100_pca.pdf")
barplot(egobp, showCategory = 20)
ggsave("BP_100_pca.pdf")

#setwd("..")
# call also for MF, BP, produce three separated plots (either use a directory for output or three separate filenames, as you prefer :))
# rbind the ego@result in a single dataframe after having addead a column 'ontology' with CC, MF or BP , use p.adjust to correct the nominal pvalue in a single run (https://www.biostars.org/p/12182/)

egoall_df <- rbind(egocc@result, egobp@result, egomf@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method='BH')
write.table(egoall_df, file = "go_results_pca_100.tsv", quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

## gsea

gsea_list <- as.data.frame(delete_genes)
gsea_list <- cbind(gsea_list, rownames(pca_genes), pca_genes$PC1)
gsea_list$`rownames(pca_genes)` <- NULL
geneList <- gsea_list$`pca_genes$PC1`
names(geneList) <- as.character(gsea_list$delete_genes)
geneList <- sort(geneList, decreasing = TRUE)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)


write.table(em@result, file = "results_GSEA_100_pca_H.tsv", quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

ridgeplot(em, showCategory = 20)
ggsave("GSEA_ridgeplot_100_pca.pdf")


##heatmap 

#delete_genes <- cbind(delete_genes, head(pca_genes$PC1, 100))
stromal <- sg %>% mutate(V1 = gsub("H_", "", V1))

vsd_f <- "/scratch/trcanmed//DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
vsd <- vsd %>% remove_rownames %>% column_to_rownames(var="genes")
vsd_f <- vsd[rownames(vsd) %in% stromal$V1,]

meda <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_lmo <- filter(meda_f, grepl("LMO_BASALE", type))
meda_lmx <- filter(meda_f, grepl("LMX_BASALE", type))
meda_lmh <- filter(meda_f, grepl("LMH", type))

meda_fin <- rbind(meda_lmx, meda_lmh)
vsd_f <- vsd_f[, colnames(vsd_f) %in% c(meda_fin$sample_id_R, "symbol")]

annotation <- as.data.frame(substr(colnames(vsd_f), 8, 10))
annotation$cases <- colnames(vsd_f)
annotation <- annotation %>% remove_rownames %>% column_to_rownames(var="cases")
colnames(annotation)[1]<-"Type"
#annotation <- annotation[order("Type")]
#t_annot <- as.data.frame(t(annotation))
#t_annot <- t_annot %>% row_to_names(t_annot, row_number = 2, remove_row = TRUE, remove_rows_above = FALSE)
t_vsd <- as.data.frame(t(vsd_f))
t_vsd$type <- as.character(substr(colnames(vsd_f), 8, 10))
t_vsd_fin <- t_vsd[order(t_vsd$type),]
t_vsd_fin$type <- NULL



#header.true <- function(df, n) {
  #names(df) <- as.character(unlist(df[n,]))
  #df[-n,]
#}
#t_annot <- header.true(t_annot, 2)




pheatmap(t_vsd, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = TRUE, show_rownames = FALSE,
         annotation_row = annotation)







