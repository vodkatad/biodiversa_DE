### get significative go padj < 0.05

lmolmxup <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE_up.tsv"
lmolmxdown <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE_down.tsv"
lmolmhup <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMH_up.tsv"
lmolmhdown <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMH_down.tsv"
lmxlmhup <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMX_BASALE.vs.LMH_up.tsv"
lmxlmhdown <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMX_BASALE.vs.LMH_down.tsv"

get_go <- function(df) {
  df <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df <- df %>% filter(p.adjust < 0.05)
  df <- df[,c(1,6)]
  rownames(df) <- NULL
  return(df)
}


lmolmxup <- get_go(lmolmxup)
lmolmxdown <- get_go(lmolmxdown)
lmolmhup <- get_go(lmolmhup)
lmolmhdown <- get_go(lmolmhdown)
lmxlmhup <- get_go(lmxlmhup)
lmxlmhdown <- get_go(lmxlmhdown)

setwd("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/")

write.table(lmolmxup, file = "go_sign_lmolmxup.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(lmolmxdown, file = "go_sign_lmolmxdown.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(lmolmhup, file = "go_sign_lmolmhup.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(lmolmhdown, file = "go_sign_lmolmhdown.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(lmxlmhup, file = "go_sign_lmxlmhup.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(lmxlmhdown, file = "go_sign_lmxlmhdown.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

