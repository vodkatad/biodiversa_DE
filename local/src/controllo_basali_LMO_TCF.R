## controllo LMO biobanca

#gene_list <- c("SPHK1","SPRR2A","NOS2","IL6","HLA-E",
#               "NLRP6","DHX58","PGLYRP4","HSPB1","CCL16",
#               "CSF2RB","LGALS9","S100A9","DDX60","IL7R","TRAV27","IGHA1",
#               "OASL","NODAL","IRF7","IL2RG","THBS1","IL4R","NEAT1","RCAN1","TMIGD1")
wt <- "/home/egrassi/onlywt_genes_go_00005.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gowt <- wt$gene

meda <- "/mnt/cold1//snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- read.table(meda, sep='\t', quote="", header=TRUE)
meda$sample_id <- gsub('-', '.', meda$sample_id_R, fixed = TRUE)
meda <- meda %>%
  filter(type %in% c("LMO_BASALE", "LMO_BASALE.1", "LMO_BASALE.2"))
meda$ssample <- substr(meda$sample_id, 1, 7)
longmodel <- meda$sample_id_R

expr <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm.tsv.gz"
expr <- read.table(expr, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
expr$genes <- gsub("H_", "", rownames(expr))
expr <- expr %>% filter(genes %in% gowt)
rownames(expr) <- expr$genes
expr$genes <- NULL
expr <- as.data.frame(t(expr))
expr$genealogy <- rownames(expr)
expr <- expr %>% filter(genealogy %in% longmodel)
expr$genealogy <- NULL

write.xlsx(expr, "espressione_geni_infiammazione_gowtprivate_lmo_basali.xlsx", rowNames=TRUE)

lfc <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv"
lfc <- read.table(lfc, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lfc$genes <- gsub("H_","", rownames(lfc))
lfc <- lfc %>% filter(genes %in% gene_list)
rownames(lfc) <- lfc$genes 
lfc$genes <- NULL

write.xlsx(lfc, "lfc_lmo_like_tcf_mut.vs.wt.xlsx", rowNames=TRUE)
