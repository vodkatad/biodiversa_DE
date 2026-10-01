## check DEG sens other

sens <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/qsens/sens_treat_cutoff0.05-TRAM.vs.NT.deseq2.tsv"
sens <- read.table(sens, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sens$genes <- rownames(sens)
sens <- sens[,c("genes", "log2FoldChange", "padj")]
colnames(sens) <- c("genes", "log2FoldChange_sens", "padj_sens")
sens_sign <- sens %>% filter(padj_sens < 0.05)
sens_sign <- sens_sign %>% filter(abs(log2FoldChange_sens) > 0.58)
sens_sign$genes <- rownames(sens_sign)
sens_sign <- sens_sign[,c("genes","log2FoldChange_sens")]

sens_up <- sens %>% filter(padj_sens < 0.05)
sens_up <- sens_up %>% filter(log2FoldChange_sens > 0.58)

sens_down <- sens %>% filter(padj_sens < 0.05)
sens_down <- sens_down %>% filter(log2FoldChange_sens < - 0.58)

other <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/qother/other_treat_cutoff0.05-TRAM.vs.NT.deseq2.tsv"
other <- read.table(other, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
other$genes <- rownames(other)
other <- other[,c("genes", "log2FoldChange", "padj")]
colnames(other) <- c("genes", "log2FoldChange_other", "padj_other")
other_sign <- other %>% filter(padj_other < 0.05)
other_sign <- other_sign %>% filter(abs(log2FoldChange_other) > 0.58)
other_sign$genes <- rownames(other_sign)
other_sign <- other_sign[,c("genes", "log2FoldChange_other")]

other_up <- other %>% filter(padj_other < 0.05)
other_up <- other_up %>% filter(log2FoldChange_other > 0.58)

other_down <- other %>% filter(padj_other < 0.05)
other_down <- other_down %>% filter(log2FoldChange_other < -0.58)

merged <- merge(sens_sign, other_sign, by="genes")

ggplot(merged, aes(log2FoldChange_sens, log2FoldChange_other))+geom_point()+geom_smooth(method = "lm")

merged$sign <- merged$log2FoldChange_sens*merged$log2FoldChange_other
merged2 <- merged %>% filter(sign < 0)

merged2 <- merge(sens, other, by="genes")

ggplot(merged2, aes(log2FoldChange_sens, log2FoldChange_other))+geom_point()+geom_smooth(method = "lm")+
  xlim(-6,6) + ylim(-6,6)

library(biomaRt)

anticorpi <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/Anticorpi_Marika_geni_da_UNI.xlsx"
anticorpi <- read.xlsx(anticorpi, sheet = 3)
genes_vec <- anticorpi$Gene_UNI |>
  na.omit() |>
  gsub("\\s*/\\s*", "/", x = _) |>
  strsplit("/") |>
  unlist() |>
  trimws() |>
  unique()
parte2 <- "/home/mferri/lista2marika.xlsx"
parte2 <- read.xlsx(parte2, sheet = 3)
parte2 <- unlist(strsplit(parte2$Gene.associato, ",\\s*"))

all <- unique(c(genes_vec, parte2))

merged2 <- merged2 %>% filter(genes %in% all)
rownames(merged2) <- merged2$genes

pheatmap(merged2[,c("log2FoldChange_sens", "log2FoldChange_other")], cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = T, fontsize_row = 4)

per_heatmap <- merged2
rownames(per_heatmap) <- per_heatmap$genes
per_heatmap$genes <- NULL

nes_cols <- c("log2FoldChange_sens", "log2FoldChange_other")
pval_cols <- c("padj_sens","padj_other")

limite <- max(abs(per_heatmap))
library(circlize)

mat_nes <- per_heatmap[, nes_cols]
mat_pval <- per_heatmap[, pval_cols]

sig_rows <- apply(mat_pval < 0.05, 1, any)

mat_nes  <- mat_nes[sig_rows, , drop = FALSE]
mat_pval <- mat_pval[sig_rows, , drop = FALSE]

library(ComplexHeatmap)
library(circlize)
mat_nes<-as.matrix(mat_nes)
limite <- max((mat_nes), na.rm = TRUE)
#colorRampPalette(c("blue","white","red"))(100)
#limite<-1
colori<-col_fun <- colorRamp2(
  c(-limite, 0, limite),
  c("#4575b4", "white", "#d73027")  
)
ord <- order(mat_nes[, "log2FoldChange_sens"], decreasing = TRUE)


ht_nes <- Heatmap(mat_nes,
                  name = "log2FC",
                  col = colori,
                  #row_order=ord,
                  cluster_columns = FALSE,
                  cluster_rows = TRUE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 4),
                  row_names_max_width = unit(12, "cm"),
                  heatmap_legend_param = list(
                    at = c(-limite, 0, limite),
                    labels = c(paste0("-", round(3, 3)), "0", paste0("+", round(3, 3))))
)
print(ht_nes)

mat_pval<-as.matrix(mat_pval)
pval_th <- 0.05
mat_pval_cat <- ifelse(mat_pval < pval_th, "significativo", "non significativo")
mat_pval_cat[is.na(mat_pval)] <- "NA" 

mat_pval_cat <- as.matrix(mat_pval_cat)

ordered_columns <- c(
  grep("_sens$", colnames(mat_pval_cat), value = TRUE),
  grep("_other$", colnames(mat_pval_cat), value = TRUE)
)

mat_pval_cat <- mat_pval_cat[, ordered_columns]

colori_pval <- c("significativo" = "goldenrod", "non significativo" = "lightgrey", "NA" = "white")


ht_pval <- Heatmap(mat_pval_cat,
                   name = "p.adjust",
                   cluster_rows = FALSE,
                   #row_order=ord,
                   cluster_columns = FALSE,
                   show_row_names = FALSE,
                   col = colori_pval,
                   rect_gp = gpar(col = "white", lwd = 1),
                   heatmap_legend_param = list(
                     at = c("significativo", "non significativo", "NA"),
                     labels = c(paste0("<", pval_th), paste0(">= ", pval_th), "NA")
                   ))

plot<-draw(ht_nes + ht_pval, heatmap_legend_side = "right")


write.table(merged2[c("genes")], file="prova.tsv", quote = FALSE, row.names = FALSE)

## espressione

fpkm <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/TRAM.vs.NT/fpkm.tsv.gz"
fpkm <- read.table(fpkm, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

calc_logFC_all <- function(df, PC = 1) {
  modelli <- unique(sub("_(NT|TRAM)_.*", "", colnames(df)))
  logFC_list <- list()
  for (m in modelli) {
    nt_cols <- grep(paste0("^", m, "_NT_"), colnames(df), value = TRUE)
    tram_cols <- grep(paste0("^", m, "_TRAM_"), colnames(df), value = TRUE)
    if(length(nt_cols) == 0 | length(tram_cols) == 0) next
    nt_mat <- df[, nt_cols, drop = FALSE]
    tram_mat <- df[, tram_cols, drop = FALSE]
    mean_nt <- rowMeans(nt_mat, na.rm = TRUE)
    mean_tram <- rowMeans(tram_mat, na.rm = TRUE)
    logFC_list[[m]] <- log2(mean_tram + PC) - log2(mean_nt + PC)
  }
  logFC_df <- do.call(cbind, logFC_list)
  rownames(logFC_df) <- rownames(df)
  return(logFC_df)
}

res <- calc_logFC_all(fpkm, PC = 1)
res <- as.data.frame(res)
res$genes <- rownames(res)
res <- res %>% filter(genes %in% genes_vec)
res$genes <- NULL

order <- "/mnt/cold1/snaketree/prj/pdxopedia/local/share/data/score_x_Ele_update_def_20260320_tram100.tsv"
order <- read.table(order, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mod_order <- order$smodel

res_ord <- res[, mod_order]

pheatmap(res_ord, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

mat_nes<-as.matrix(res_ord)
limite <- max((mat_nes), na.rm = TRUE)
#colorRampPalette(c("blue","white","red"))(100)
#limite<-1
colori<-col_fun <- colorRamp2(
  c(-4, 0, 4),
  c("#4575b4", "white", "#d73027")  
)

ht_nes <- Heatmap(mat_nes,
                  col = colori,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(12, "cm"),
                  heatmap_legend_param = list(
                    at = c(-4, 0, 4),
                    labels = c(paste0("-", round(4, 4)), "0", paste0("+", round(4, 4))))
)
print(ht_nes)

make_mean_matrix <- function(df) {
  modelli <- unique(sub("_(NT|TRAM)_.*", "", colnames(df)))
  res_list <- list()
  for (m in modelli) {
    nt_cols <- grep(paste0("^", m, "_NT_"), colnames(df), value = TRUE)
    tram_cols <- grep(paste0("^", m, "_TRAM_"), colnames(df), value = TRUE)
    # NT
    if (length(nt_cols) > 0) {
      nt_mat <- df[, nt_cols, drop = FALSE]
      nt_mean <- rowMeans(nt_mat, na.rm = TRUE)
      res_list[[paste0(m, "_NT")]] <- nt_mean
    }
    # TRAM
    if (length(tram_cols) > 0) {
      tram_mat <- df[, tram_cols, drop = FALSE]
      tram_mean <- rowMeans(tram_mat, na.rm = TRUE)
      res_list[[paste0(m, "_TRAM")]] <- tram_mean
    }
  }
  res <- do.call(cbind, res_list)
  rownames(res) <- rownames(df)
  return(res)
}

mean_mat <- make_mean_matrix(fpkm)
mean_mat <- as.data.frame(mean_mat)
log_mat <- log2(mean_mat + 1)
log_mat$genes <- rownames(log_mat)
log_mat <- log_mat %>% filter(genes %in% genes_vec)
log_mat$genes <- NULL

pheatmap(log_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

col_order <- unlist(lapply(mod_order, function(m) {
  c(paste0(m, "_NT"), paste0(m, "_TRAM"))
}))

log_mat_ord <- log_mat[, col_order]

pheatmap(log_mat_ord, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

log_mat_scaled <- t(scale(t(log_mat_ord)))

pheatmap(log_mat_scaled, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

mat_nes<-as.matrix(log_mat_scaled)
limite <- max((mat_nes), na.rm = TRUE)
#colorRampPalette(c("blue","white","red"))(100)
#limite<-1
colori<-col_fun <- colorRamp2(
  c(-7, 0, 7),
  c("#4575b4", "white", "#d73027")  
)

ht_nes <- Heatmap(mat_nes,
                  col = colori,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(12, "cm"),
                  heatmap_legend_param = list(
                    at = c(-7, 0, 7),
                    labels = c(paste0("-", round(7, 7)), "0", paste0("+", round(7, 7))))
)
print(ht_nes)

## check NT per TGFB2 e CDKN2B

expr <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/TRAM.vs.NT/fpkm.tsv.gz"
expr <- read.table(expr, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
expr_nt <- expr[, grepl("_NT_", colnames(expr))]
patient_ids <- sub("^(CRC[0-9]+).*", "\\1", colnames(expr_nt))
expr_nt_mean <- sapply(unique(patient_ids), function(id) {
  rowMeans(expr_nt[, patient_ids == id, drop = FALSE], na.rm = TRUE)
})

expr_nt_mean <- as.data.frame(expr_nt_mean)
expr_nt_mean$genes <- rownames(expr_nt_mean)

#find <- c("TGFB2", "CDKN2B")
find <- c("MYC",
          "CCND1",
          "CDK4",
          "CDK6",
          "CDK2",
          "CCNE1",
          "CDKN2A",
          "CDKN2B",
          "CDKN1A",
          "CDKN1B",
          "CDKN2C",
          "CDKN2D",
          "TP53",
          "APC",
          "PTEN",
          "SMAD2",
          "SMAD3",
          "SMAD4",
          "TGFBR1",
          "TGFBR2",
          "CTNNB1",
          "TGFB1",
          "TGFB2",
          "ACVR1B",
          "ACVR2A")

expr_nt_mean <- expr_nt_mean %>% filter(genes %in% find)
expr_nt_mean$genes <- NULL
expr_nt_mean_log2 <- log2(expr_nt_mean + 1)
expr_nt_mean_log2_order <- expr_nt_mean_log2[, mod_order]

mat_nes<-as.matrix(expr_nt_mean_log2_order)
limite <- max((mat_nes), na.rm = TRUE)
#colorRampPalette(c("blue","white","red"))(100)
#limite<-1
colori<-col_fun <- colorRamp2(
  c(-10, 0, 10),
  c("#4575b4", "white", "#d73027")  
)

ht_nes <- Heatmap(mat_nes,
                  col = colori,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(12, "cm"),
                  heatmap_legend_param = list(
                    at = c(-10, 0, 10),
                    labels = c(paste0("-", round(10, 10)), "0", paste0("+", round(10, 10))))
)
print(ht_nes)

write.xlsx(expr_nt_mean_log2_order, file="/home/mferri/expr_geni_Marika_mean_log2_order.xlsx", rowNames=T, colNames=T)

## wilcoxon per geni anticorpi coi casi di sens e other

fpkm <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/TRAM.vs.NT/fpkm.tsv.gz"
fpkm <- read.table(fpkm, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

calc_logFC_all <- function(df, PC = 1) {
  modelli <- unique(sub("_(NT|TRAM)_.*", "", colnames(df)))
  logFC_list <- list()
  for (m in modelli) {
    nt_cols <- grep(paste0("^", m, "_NT_"), colnames(df), value = TRUE)
    tram_cols <- grep(paste0("^", m, "_TRAM_"), colnames(df), value = TRUE)
    if(length(nt_cols) == 0 | length(tram_cols) == 0) next
    nt_mat <- df[, nt_cols, drop = FALSE]
    tram_mat <- df[, tram_cols, drop = FALSE]
    mean_nt <- rowMeans(nt_mat, na.rm = TRUE)
    mean_tram <- rowMeans(tram_mat, na.rm = TRUE)
    logFC_list[[m]] <- log2(mean_tram + PC) - log2(mean_nt + PC)
  }
  logFC_df <- do.call(cbind, logFC_list)
  rownames(logFC_df) <- rownames(df)
  return(logFC_df)
}

res <- calc_logFC_all(fpkm, PC = 1)
res <- as.data.frame(res)
res$genes <- rownames(res)
res <- res %>% filter(genes %in% genes_vec)
res$genes <- NULL

## prendere solo i casi che ci interessano

casisens <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/qsens/sens_samples_data"
casisens <- read.table(casisens, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casisens <- unique(casisens$smodel)

casiother <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/qother/other_samples_data"
casiother <- read.table(casiother, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casiother <- unique(casiother$smodel)

casiperwilcoxon <- c(casiother, casisens)

res_sub <- res[, colnames(res) %in% casiperwilcoxon, drop = FALSE]

sens_models  <- intersect(casisens, colnames(res_sub))
other_models <- intersect(casiother, colnames(res_sub))

results <- apply(res_sub, 1, function(x) {
  wilcox.test(
    x[sens_models],
    x[other_models],
    alternative = "two.sided"
  )$p.value
})

results <- data.frame(
  gene = rownames(res_sub),
  pval = results
)

results <- results[order(results$pval), ]
results$padj <- p.adjust(results$pval, method = "BH")

## heatmap mutazione per geni su cui è stata fatta l'analisi GISTIC

mut <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/Marika/Marika_check_mut_del/Marika_check/Marika_check_specific-mutations_list_GOI_PDO_alltiers.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

genes <- unique(mut$mutated_gene)

full_grid <- expand.grid(
  mutated_gene = genes,
  smodel = casiperwilcoxon
)

mut_bin <- mut %>%
  distinct(smodel, mutated_gene) %>%
  mutate(value = 1)

heat_df <- full_grid %>%
  left_join(mut_bin,
            by = c("mutated_gene", "smodel")) %>%
  mutate(value = ifelse(is.na(value), 0, value))

mat <- heat_df %>%
  pivot_wider(
    names_from = smodel,
    values_from = value
  )

mat <- as.data.frame(mat)

rownames(mat) <- mat$mutated_gene
mat$mutated_gene <- NULL

mat <- as.matrix(mat)

ordered_samples <- c(sens_models, other_models)

mat <- mat[, ordered_samples]

group_annotation <- data.frame(
  deg = c(
    rep("sens", length(sens_models)),
    rep("other", length(other_models))
  )
)

rownames(group_annotation) <- ordered_samples

ha <- HeatmapAnnotation(
  deg = group_annotation$deg,
  col = list(
    deg = c(
      sens = "#f4a09c",
      other = "#21d5df"
    )
  )
)

Heatmap(
  mat,
  name = "Mutation",
  col = c(
    "0" = "grey80",
    "1" = "#ff6a4d"
  ),
  top_annotation = ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  rect_gp = gpar(col = "grey70"),
  row_names_side = "right"
)
