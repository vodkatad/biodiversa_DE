### first part analysis keratine for AA deg

setwd("/scratch/trcanmed/DE_RNASeq/dataset/connector_Ac_all/")

### FARE TABLE PER KRT E EMT

geni_krt <- "../connector_for_all/keratine_padre.tsv"
geni_krt <- read.table(geni_krt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
krt <- geni_krt$SYMBOL

aa <- "Ac_Aa/col_cutoff0.05-Ac.vs.Aa.deseq2_Ac_Aa.tsv"
ab <- "Ac_Ab/col_cutoff0.05-Ac.vs.Ab.deseq2_Ac_Ab.tsv"
ba <- "Ac_Ba/col_cutoff0.05-Ac.vs.Ba.deseq2_Ac_Ba.tsv"
bb <- "Ac_Bb/col_cutoff0.05-Ac.vs.Bb.deseq2_Ac_Bb.tsv"
bc <- "Ac_Bc/col_cutoff0.05-Ac.vs.Bc.deseq2_Ac_Bc.tsv"
c <- "Ac_C/col_cutoff0.05-Ac.vs.C.deseq2_Ac_C.tsv"
#pd_or <- "/scratch/trcanmed/DE_RNASeq/dataset/recist_connector_DEG/classification_cutoff0.05-PD.vs.OR.deseq2.tsv"

get_file_ready <- function(df, what_class) {
  df_f <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  genes_h <- rownames(df_f)
  rownames(df_f) <- NULL
  df_f <- cbind(genes_h,df_f)
  df_f <- df_f %>% mutate(genes_h = gsub("H_", "", genes_h))
  colnames(df_f)[colnames(df_f) == 'genes_h'] <- 'symbol'
  df_f <- df_f %>% filter(symbol %in% what_class)
  return(df_f)
}

aak <- get_file_ready(aa, krt)
abk <- get_file_ready(ab, krt)
bak <- get_file_ready(ba, krt)
bbk <- get_file_ready(bb, krt)
bck <- get_file_ready(bc, krt)
ck <- get_file_ready(c, krt)

get_column_mean <- function(df) {
  df_mean <- mean(df$log2FoldChange)
  return(df_mean)
}

mean_aa <- get_column_mean(aak)
mean_ab <- get_column_mean(abk)
mean_ba <- get_column_mean(bak)
mean_bb <- get_column_mean(bbk)
mean_bc <- get_column_mean(bck)
mean_c <- get_column_mean(ck)
#mean_pd_or <- -(get_column_mean(pd_or))

means_tot <- as.data.frame(cbind(c("Aa", "Ab", "Ba", "Bb", "Bc", "C")))
means_tot <- cbind(means_tot, c(mean_aa, mean_ab, mean_ba, mean_bb, mean_bc, mean_c))
colnames(means_tot) <- c("Ac_vs", "mean_LFC_GO_sign_keratine_large")

ggplot(data=means_tot, aes(x=Ac_vs, y=mean_LFC_GO_sign_keratine_large)) + geom_bar(stat="identity")

geni_emt <- "../connector_for_all/epi_to_mes.tsv"
geni_emt <- read.table(geni_emt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
emt <- geni_emt$SYMBOL

aae <- get_file_ready(aa, emt)
abe <- get_file_ready(ab, emt)
bae <- get_file_ready(ba, emt)
bbe <- get_file_ready(bb, emt)
bce <- get_file_ready(bc, emt)
ce <- get_file_ready(c, emt)

get_column_mean <- function(df) {
  df_mean <- mean(df$log2FoldChange)
  return(df_mean)
}

mean_aae <- get_column_mean(aae)
mean_abe <- get_column_mean(abe)
mean_bae <- get_column_mean(bae)
mean_bbe <- get_column_mean(bbe)
mean_bce <- get_column_mean(bce)
mean_ce <- get_column_mean(ce)
#mean_pd_or <- -(get_column_mean(pd_or))

means_tote <- as.data.frame(cbind(c("Aa", "Ab", "Ba", "Bb", "Bc", "C")))
means_tote <- cbind(means_tote, c(mean_aae, mean_abe, mean_bae, mean_bbe, mean_bce, mean_ce))
colnames(means_tote) <- c("Ac_vs", "mean_LFC_GO_sign_emt_small")

ggplot(data=means_tote, aes(x=Ac_vs, y=mean_LFC_GO_sign_emt_small)) + geom_bar(stat="identity")

geni_emtl <- "../connector_for_all/emt_large.tsv"
geni_emtl <- read.table(geni_emtl, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
emtl <- geni_emtl$SYMBOL

aael <- get_file_ready(aa, emtl)
abel <- get_file_ready(ab, emtl)
bael <- get_file_ready(ba, emtl)
bbel <- get_file_ready(bb, emtl)
bcel <- get_file_ready(bc, emtl)
cel <- get_file_ready(c, emtl)

get_column_mean <- function(df) {
  df_mean <- mean(df$log2FoldChange)
  return(df_mean)
}

mean_aael <- get_column_mean(aael)
mean_abel <- get_column_mean(abel)
mean_bael <- get_column_mean(bael)
mean_bbel <- get_column_mean(bbel)
mean_bcel <- get_column_mean(bcel)
mean_cel <- get_column_mean(cel)
#mean_pd_or <- -(get_column_mean(pd_or))

means_totel <- as.data.frame(cbind(c("Aa", "Ab", "Ba", "Bb", "Bc", "C")))
means_totel <- cbind(means_totel, c(mean_aael, mean_abel, mean_bael, mean_bbel, mean_bcel, mean_cel))
colnames(means_totel) <- c("Ac_vs", "mean_LFC_GO_sign_emt_large")

ggplot(data=means_totel, aes(x=Ac_vs, y=mean_LFC_GO_sign_emt_large)) + geom_bar(stat="identity")

aau <- "Ac_Aa/GO_results_col_cutoff0.05-Ac.vs.Aa_up_Ac_Aa.tsv"
aad <- "Ac_Aa/GO_results_col_cutoff0.05-Ac.vs.Aa_down_Ac_Aa.tsv"
abu <- "Ac_Ab/GO_results_col_cutoff0.05-Ac.vs.Ab_up_Ac_Ab.tsv"
abd <- "Ac_Ab/GO_results_col_cutoff0.05-Ac.vs.Ab_down_Ac_Ab.tsv"
bau <- "Ac_Ba/GO_results_col_cutoff0.05-Ac.vs.Ba_up_Ac_Ba.tsv"
bad <- "Ac_Ba/GO_results_col_cutoff0.05-Ac.vs.Ba_down_Ac_Ba.tsv"
bbu <- "Ac_Bb/GO_results_col_cutoff0.05-Ac.vs.Bb_up_Ac_Bb.tsv"
bbd <- "Ac_Bb/GO_results_col_cutoff0.05-Ac.vs.Bb_down_Ac_Bb.tsv"
bcu <- "Ac_Bc/GO_results_col_cutoff0.05-Ac.vs.Bc_up_Ac_Bc.tsv"
bcd <- "Ac_Bc/GO_results_col_cutoff0.05-Ac.vs.Bc_down_Ac_Bc.tsv"
cu <- "Ac_C/GO_results_col_cutoff0.05-Ac.vs.C_up_Ac_C.tsv"
cd <- "Ac_C/GO_results_col_cutoff0.05-Ac.vs.C_down_Ac_C.tsv"

get_go_sign <- function(df) {
  df_f <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df_f <- df_f %>% filter(p.adjust < 0.05)
  return(df_f)
}

goacu <- get_go_sign(aau)
goacd <- get_go_sign(aad)
#goabu <- get_go_sign(abu)
#goabd <- get_go_sign(abd)
gobau <- get_go_sign(bau)
gobad <- get_go_sign(bad) #sign GO_Ac_Ba_down
gobbu <- get_go_sign(bbu) 
gobbd <- get_go_sign(bbd) #sign Go_Ac_Bb_down
gobcu <- get_go_sign(bcu)
gobcd <- get_go_sign(bcd)
gocu <- get_go_sign(cu)
gocd <- get_go_sign(cd)

saa <- "Ac_Aa/samples_datas_Ac_Aa"
sab <- "Ac_Ab/samples_datas_Ac_Ab"
sba <- "Ac_Ba/samples_datas_Ac_Ba"
sbb <- "Ac_Bb/samples_datas_Ac_Bb"
sbc <- "Ac_Bc/samples_datas_Ac_Bc"
sc <- "Ac_C/samples_datas_Ac_C"

tell_me_table <- function(df) {
  df_f <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  return(table(df_f$col))
}

saa <- tell_me_table(saa)
sab <- tell_me_table(sab)
sba <- tell_me_table(sba)
sbb <- tell_me_table(sbb)
sbc <- tell_me_table(sbc)
sc <- tell_me_table(sc)
