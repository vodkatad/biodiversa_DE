## obtain lfc for all genes 

## make pseudocounts for all singles deg

directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/fpkm_singoli_N2.vs.NE_def"

tsv_files <- list.files(directory_path, pattern = "\\.gz$", full.names = TRUE)

data_list <- list()
for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  data_list[[file_name]] <- df
}

colnames(data_list[["CRC0542_fpkm.tsv"]]) <- c("CRC0542_NE_R1", "CRC0542_N2_R1", "CRC0542_NE_R2", "CRC0542_N2_R2")

calc_logFC_df <- function(df, sample_name, PC = 1) { #PC = 0.0001
  cols_N2 <- grep("N2", colnames(df), value = TRUE)
  cols_NE <- grep("NE", colnames(df), value = TRUE)
  if(length(cols_N2) != length(cols_NE)) {
    stop("Numero di colonne N2 e NE diverso!")
  }
  logFC_matrix <- mapply(function(n2, ne) {
    log2((df[[n2]] + PC) / (df[[ne]] + PC))
  }, cols_N2, cols_NE)
  logFC <- rowMeans(logFC_matrix)
  df_logFC <- data.frame(row.names = rownames(df))
  df_logFC[[paste0(sample_name, "_LFC")]] <- logFC
  return(df_logFC)
}

logFC_list <- mapply(calc_logFC_df, data_list, names(data_list), SIMPLIFY = FALSE)
names(logFC_list) <- gsub("_fpkm\\.tsv$", "", names(logFC_list))


all_row_names <- unique(unlist(lapply(logFC_list, rownames)))
template_df <- data.frame(row.names = all_row_names)
for (name in names(logFC_list)) {
  df <- logFC_list[[name]]
  template_df[[paste0(name, "_LFC")]] <- df[match(all_row_names, rownames(df)), , drop = TRUE]
}

res <- template_df
res$genes <- rownames(res)

lfc_res <- res

## parte heatmap

res <- lfc_res
bcat_all <- "/home/mferri/GO0016055_wnt_signaling_pathway.tsv"
bcat_all <- read.table(bcat_all, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bcat_all_genes <- unique(bcat_all$SYMBOL)

res <- res %>% filter(genes %in% bcat_all_genes)
res$genes <- NULL
colnames(res) <- gsub("_LFC", "", colnames(res))

## annotatio col 
order_dep_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx" 
order_dep <- read_xlsx(order_dep_f) 
order_dep$quartile <- ntile(order_dep$co_comp_N2, 4) 
names(order_dep)[names(order_dep)=="case"] <- "model" 
order_dep <- as.data.frame(order_dep) 
rownames(order_dep) <- order_dep$model 
mut <- c("CRC0148", "CRC1331", "CRC0399", "CRC1278", "CRC0277", 
         "CRC0327", "CRC1729", "CRC0152", "CRC0196", "CRC0059", 
         "CRC0065", "CRC0464", "CRC0316", "CRC1239") 
for (i in rownames(order_dep)) { 
  if (order_dep[i, "model"] %in% mut) { order_dep[i, "mut"] <- "MUT" 
  } else { order_dep[i, "mut"] <- "WT" 
  } 
} 
order_dep$co_comp_N2 <- NULL 
order_dep <- order_dep %>% filter(!model == c("CRC0743", "CRC0152")) 
order_dep$model <- NULL 
order_dep <- order_dep[order(order_dep$mut),] 
#order_dep$quartile <- NULL 
order_depnoq <- order_dep 
order_depnoq$quartile <- NULL

res_ordered  <- res[, rownames(order_dep), drop = FALSE]

d <- as.matrix(res_ordered)

minv <- -1.5
maxv <- 1.5
neutral_value <- 0
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue","lightblue"))(n = length(bk1)-1),
                "#e1e1e1","#e1e1e1",
                colorRampPalette(colors = c("tomato1","darkred"))(n = length(bk2)-1))


pheatmap(d, cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = bk, color = my_palette,
         annotation_col = order_dep, na_col = "#FFFFFF")


# scatter <- as.data.frame(t(scatter))
# scatter$row_mean <- rowMeans(scatter, na.rm = TRUE)
# order_dep$model <- rownames(order_dep)
# scatter$model <- rownames(scatter)
# scatter <- scatter[,c("model", "row_mean")]
# scatter <- merge(scatter, order_dep, by="model")
# 
# ggplot(scatter, aes(x=row_mean, y=co_comp_N2, color=mut))+geom_point()+geom_smooth(method = lm)+geom_text_repel(aes(label = model))+xlab("Mean_LFC_DEG_model")

## segnare la signficatività a lato della heatmap 
## mut a sinistra e wt a destra

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut$genes <- rownames(mut)
mut <- mut %>% filter(genes %in% bcat_all_genes)
for (i in rownames(mut)) {
  if (mut[i, "padj"] < 0.05) {
    mut[i, "sign"] <- "sign"
  } else {
    mut[i, "sign"] <- "no"
  }
}
mut <- mut[, c("genes", "sign")]

wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE//WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt$genes <- rownames(wt)
wt <- wt %>% filter(genes %in% bcat_all_genes)
for (i in rownames(wt)) {
  if (wt[i, "padj"] < 0.05) {
    wt[i, "sign"] <- "sign"
  } else {
    wt[i, "sign"] <- "no"
  }
}

wt <- wt[, c("genes", "sign")]

# # geni che mancano nei deg singoli
# new_genes <- c("HESX1", "DKKL1")
# 
# new_rows <- matrix(NA, nrow = length(new_genes), ncol = ncol(d))
# rownames(new_rows) <- new_genes
# colnames(new_rows) <- colnames(d)
# 
# d <- rbind(d, new_rows)

col_fun <- colorRamp2(c(-1, 0, 1),
                      c("darkblue", "#e1e1e1", "darkred"))

row_anno_mut <- mut[match(rownames(d), mut$genes), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

top_anno <- HeatmapAnnotation(df = order_dep)

ht <- Heatmap(
  as.matrix(d),
  name = "LFC",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  na_col = "#FFFFFF",
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(MUT = row_anno_mut,
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)

draw(ht, merge_legend = TRUE)

geni_selected <- rownames(d)
## add new annotation
bcat_all$GO <- "WNT_signaling_pathway"
bcat_all <- bcat_all[,c("SYMBOL", "GO")]
bcat_all <- bcat_all %>% filter(SYMBOL %in% geni_selected)
bcat_all <- bcat_all[!duplicated(bcat_all$SYMBOL),]

positivi_wtn <- "/home/mferri/GO0030177_positive_regulation_wnt.tsv"
wntplus <- read.table(positivi_wtn, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wntplus$GO <- "WNT_positive_regulation"
wntplus <- wntplus[!duplicated(wntplus$SYMBOL),]
wntplus <- wntplus[,c("SYMBOL", "GO")]

WNT <- merge(bcat_all, wntplus, by="SYMBOL", all.x = TRUE)

negative_wnt <- "/home/mferri/GO0030178_negative_regulation_wnt.tsv"
wntm <- read.table(negative_wnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wntm$GO <- "WNT_negative_regulation"
wntm <- wntm[!duplicated(wntm$SYMBOL),]
wntm <- wntm[,c("SYMBOL", "GO")]

WNT <- merge(WNT, wntm, by="SYMBOL", all.x=TRUE)

canonical_wnt <- "/home/mferri/GO0060070_canonical_wnt_signaling_pathway.tsv"
cwnt <- read.table(canonical_wnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cwnt$GO <- "Canonical_Wnt_signaling_pathway"
cwnt <- cwnt[!duplicated(cwnt$SYMBOL),]
cwnt <- cwnt[,c("SYMBOL", "GO")]

WNT <- merge(WNT, cwnt, by="SYMBOL", all.x=TRUE)

rownames(WNT) <- WNT$SYMBOL
WNT$SYMBOL <- NULL
colnames(WNT) <- c("WNT_signaling_pathway", "WNT_positive_regulation", "WNT_negative_regulation", "Canonical_Wnt_signaling_pathway")
WNT <- WNT[,c(1,4,2,3)]

## new heatmap
wnt_df <- WNT[match(rownames(d), rownames(WNT)), ]

col_wnt <- list(
  WNT_signaling_pathway = c("WNT_signaling_pathway" = "#d73027", "NA" = "grey90"),
  Canonical_Wnt_signaling_pathway = c("Canonical_Wnt_signaling_pathway" = "#4575b4", "NA" = "grey90"),
  WNT_positive_regulation = c("WNT_positive_regulation" = "#1a9850", "NA" = "grey90"),
  WNT_negative_regulation = c("WNT_negative_regulation" = "#fee08b", "NA" = "grey90")
)

left_anno <- rowAnnotation(
  MUT = row_anno_mut,
  WT  = row_anno_wt,
  WNT_signaling_pathway = wnt_df$WNT_signaling_pathway,
  Canonical_Wnt_signaling_pathway = wnt_df$Canonical_Wnt_signaling_pathway,
  WNT_positive_regulation = wnt_df$WNT_positive_regulation,
  WNT_negative_regulation = wnt_df$WNT_negative_regulation,
  col = c(
    list(
      MUT = c("sign" = "tomato", "no" = "grey80"),
      WT  = c("sign" = "steelblue", "no" = "grey80")
    ),
    col_wnt
  ),
  annotation_name_side = "top"
)

ht <- Heatmap(
  as.matrix(d),
  name = "LFC",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = col_fun,
  na_col = "#FFFFFF",
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 3),
  top_annotation = top_anno,
  left_annotation = left_anno
)

draw(ht, merge_legend = TRUE)

