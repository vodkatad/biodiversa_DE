### first part analysis keratine for AA deg

setwd("/scratch/trcanmed/DE_RNASeq/dataset/connector_for_all/")

geni <- "QuickGO-annotations-1655821430897-20220621.tsv"
geni <- read.table(geni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- geni$SYMBOL

ab <- "Aa_Ab/col_cutoff0.05-Aa.vs.Ab.deseq2_Aa_Ab.tsv"
ac <- "Aa_Ac/col_cutoff0.05-Aa.vs.Ac.deseq2.tsv"
ba <- "Aa_Ba/col_cutoff0.05-Aa.vs.Ba.deseq2_Aa_Ba.tsv"
bb <- "Aa_Bb/col_cutoff0.05-Aa.vs.Bb.deseq2_Aa_Bb.tsv"
bc <- "Aa_Bc/col_cutoff0.05-Aa.vs.Bc.deseq2_Aa_Bc.tsv"
c <- "Aa_C/col_cutoff0.05-Aa.vs.C.deseq2_Aa_C.tsv"
pd_or <- "/scratch/trcanmed/DE_RNASeq/dataset/recist_connector_DEG/classification_cutoff0.05-PD.vs.OR.deseq2.tsv"

get_column_mean <- function(df) {
  df_f <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  genes_h <- rownames(df_f)
  rownames(df_f) <- NULL
  df_f <- cbind(genes_h,df_f)
  df_f <- df_f %>% mutate(genes_h = gsub("H_", "", genes_h))
  colnames(df_f)[colnames(df_f) == 'genes_h'] <- 'symbol'
  df_f <- df_f %>% filter(symbol %in% genes)
  df_mean <- mean(df_f$log2FoldChange)
  return(df_mean)
}

mean_ab <- get_column_mean(ab)
mean_ac <- get_column_mean(ac)
mean_ba <- get_column_mean(ba)
mean_bb <- get_column_mean(bb)
mean_bc <- get_column_mean(bc)
mean_c <- get_column_mean(c)
mean_pd_or <- -(get_column_mean(pd_or))

means_tot <- as.data.frame(cbind(c("Ab", "Ac", "Ba", "Bb", "Bc", "C", "OR_PD")))
means_tot <- cbind(means_tot, c(mean_ab, mean_ac, mean_ba, mean_bb, mean_bc, mean_c, mean_pd_or))
colnames(means_tot) <- c("Aa_vs", "mean_LFC_GO_sign")

ggplot(data=means_tot, aes(x=Aa_vs, y=mean_LFC_GO_sign)) + geom_bar(stat="identity")

abu <- "Aa_Ab/GO_results_col_cutoff0.05-Aa.vs.Ab_up_Aa_Ab.tsv"
abd <- "Aa_Ab/GO_results_col_cutoff0.05-Aa.vs.Ab_down_Aa_Ab.tsv"
acu <- "Aa_Ac/GO_results_col_cutoff0.05-Aa.vs.Ac_up.tsv"
acd <- "Aa_Ac/GO_results_col_cutoff0.05-Aa.vs.Ac_down.tsv"
bau <- "Aa_Ba/GO_results_col_cutoff0.05-Aa.vs.Ba_up_Aa_Ba.tsv"
bad <- "Aa_Ba/GO_results_col_cutoff0.05-Aa.vs.Ba_down_Aa_Ba.tsv"
bbu <- "Aa_Bb/GO_results_col_cutoff0.05-Aa.vs.Bb_up_Aa_Bb.tsv"
bbd <- "Aa_Bb/GO_results_col_cutoff0.05-Aa.vs.Bb_down_Aa_Bb.tsv"
bcu <- "Aa_Bc/GO_results_col_cutoff0.05-Aa.vs.Bc_up_Aa_Bc.tsv"
bcd <- "Aa_Bc/GO_results_col_cutoff0.05-Aa.vs.Bc_down_Aa_Bc.tsv"
cu <- "Aa_C/GO_results_col_cutoff0.05-Aa.vs.C_up_Aa_C.tsv"
cd <- "Aa_C/GO_results_col_cutoff0.05-Aa.vs.C_down_Aa_C.tsv"

get_go_sign <- function(df) {
  df_f <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df_f <- df_f %>% filter(p.adjust < 0.05)
  return(df_f)
}

goabu <- get_go_sign(abu)
goabd <- get_go_sign(abd)
goacu <- get_go_sign(acu)
goacd <- get_go_sign(acd)
gobau <- get_go_sign(bau)
gobad <- get_go_sign(bad) #sign
gobbu <- get_go_sign(bbu) 
gobbd <- get_go_sign(bbd) #sign
gobcu <- get_go_sign(bcu)
gobcd <- get_go_sign(bcd)
gocu <- get_go_sign(cu)
gocd <- get_go_sign(cd)
  
