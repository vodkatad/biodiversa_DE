## confronto GO con GSEA C5

### TCF main

go <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv"
go <- read.table(go, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go <- go %>% filter(p.adjust < 0.05)
go <- go[order(go$p.adjust),]

gsea <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GSEA_results_C5_geno_cutoff0.05-N2.vs.NE.tsv"
gsea <- read.table(gsea, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea <- gsea %>% filter(NES < 0)
gsea <- gsea %>% filter(p.adjust < 0.05)
gsea <- gsea[!grepl("HP_", gsea$ID), ]
gsea <- gsea[order(gsea$p.adjust),]

gsea_bp <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GSEA_results_C5_BP_geno_cutoff0.05-N2.vs.NE.tsv"
gsea_bp <- read.table(gsea_bp, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_bp <- gsea_bp %>% filter(NES < 0)
gsea_bp <- gsea_bp %>% filter(p.adjust < 0.05)

## TCF basali

go_b <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/GO_results_geno_mut_cutoff0.05-NE_MUT.vs.NE_WT_down.tsv"
go_b <- read.table(go_b, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_b <- go_b %>% filter(p.adjust < 0.05)

gsea_b <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/geno_mut_GSEA_results_C5_cutoff0.05-NE_MUT.vs.NE_WT.tsv"
gsea_b <- read.table(gsea_b, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_b <- gsea_b %>% filter(p.adjust < 0.05)


vsd <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/vsd.tsv.gz", quote = "",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$genes <- rownames(vsd)
vsd <- vsd %>% filter(genes == "EPHA2")
vsd$genes <- NULL
vsd <- as.data.frame(t(vsd))


