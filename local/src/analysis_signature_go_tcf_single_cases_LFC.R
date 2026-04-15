library(ComplexHeatmap)
library(circlize)

## check lfc signatures in single degs 

bcat_go <- c("ZNRF3", "KREMEN1", "NKD1", "NOTUM", "HHEX", "PTK7", "LGR5", "LGR6", "APOE",  
            "NID1", "ADGRA2", "WNT8B", "DKK3", "DLX5", "HESX1", "DACT3", "FZD10", "TERT",  
            "BARX1", "FZD9", "RSPO4", "TMEM131L", "DISC1", "WNT10B", "FZD2", "GPRC5B", "LRP4",  
            "TPBGL", "DKKL1", "TMEM88B", "ARHGEF19")
## geni bcat_progress
bcat_progress <- c("LGR5", "EPHB3", "NKD1", "ZNRF3", "ASCL2", "GINS2",
                   "MYC", "SP5", "GINS3", "CFCA4", "TEAD4", "BMP4", "RNF43", 
                   "AXIN2", "CCND1", "DKK1")

directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/logfc_singoli_N2.vs.NE_def/"

tsv_files <- list.files(directory_path, pattern = "\\.tsv$", full.names = TRUE)
data_frames <- list()
padj_frames <- list()

for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  df_lfc  <- data[, 2, drop = FALSE]
  df_padj <- data[, 6, drop = FALSE]
  colnames(df_lfc)  <- paste0(file_name, "_LFC")
  colnames(df_padj) <- paste0(file_name, "_padj")
  df_lfc$genes  <- rownames(data)
  df_padj$genes <- rownames(data)
  df_lfc  <- df_lfc  %>% filter(genes %in% bcat_progress)
  df_padj <- df_padj %>% filter(genes %in% bcat_progress)
  rownames(df_lfc)  <- df_lfc$genes;  df_lfc$genes  <- NULL
  rownames(df_padj) <- df_padj$genes; df_padj$genes <- NULL
  data_frames[[file_name]] <- df_lfc
  padj_frames[[file_name]] <- df_padj
}

all_row_names <- unique(unlist(lapply(data_frames, rownames)))
template_df   <- data.frame(row.names = all_row_names)
template_padj <- data.frame(row.names = all_row_names)

for (name in names(data_frames)) {
  df <- data_frames[[name]]
  template_df[[name]] <- df[match(all_row_names, rownames(df)), ]
  
  dfp <- padj_frames[[name]]
  template_padj[[name]] <- dfp[match(all_row_names, rownames(dfp)), ]
}

res <- template_df
res_padj <- template_padj

## togliere quei geni che non hanno lfc per almeno 15/30 campioni 
#res_filtered <- res[rowSums(!is.na(res)) >= 15, ] 
#res_filtered <- as.data.frame(t(res_filtered))
res_filtered <- res
res_padj_filtered <- res_padj[rownames(res_filtered), , drop=FALSE]

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

res_ordered  <- res_filtered[, rownames(order_dep), drop = FALSE]
padj_ordered <- res_padj_filtered[, rownames(order_dep), drop = FALSE]

d <- as.matrix(res_ordered)
padj_mat <- as.matrix(padj_ordered)

stars <- ifelse(padj_mat < 0.05, "*", "")

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
         annotation_col = order_dep, na_col = "#FFFFFF",
         display_numbers = stars, number_color = "black")


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
mut <- mut %>% filter(genes %in% bcat_progress)
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
wt <- wt %>% filter(genes %in% bcat_progress)
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

col_fun <- colorRamp2(c(-1.5, 0, 1.5),
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

analysis <- d
analysis <- as.data.frame(t(analysis))
analysis$sample <- rownames(analysis)
order_dep$sample <- rownames(order_dep)
analysis <- merge(analysis, order_dep, by="sample")
analysis$quartile <- NULL
wt <- analysis %>% filter(mut == "WT")
wt$mut <- NULL
rownames(wt) <- wt$sample
wt$sample <- NULL
wt <- as.data.frame(t(wt))
mut <- analysis %>% filter(mut == "MUT")
mut$mut <- NULL
rownames(mut) <- mut$sample
mut$sample <- NULL
mut <- as.data.frame(t(mut))

#mut <- mut[rowSums(is.na(mut)) != ncol(mut), ]
#wt <- wt[!(rownames(wt) %in% "TMEM88B"), ]
mut <- mut[rownames(wt), ]

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_bcat_progress <- data.frame(
  gene = names(results),
  value = unlist(results)
)


## jak-stat
jak <- c("PTK2B", "IL23R", "IL6R", "CCL5", "IL7R", "TNF", "MIR221", "GHR", "IFNL1", "IL6") 

directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/logfc_singoli_N2.vs.NE_def/"

tsv_files <- list.files(directory_path, pattern = "\\.tsv$", full.names = TRUE)
data_frames <- list()
padj_frames <- list()

for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  df_lfc  <- data[, 2, drop = FALSE]
  df_padj <- data[, 6, drop = FALSE]
  colnames(df_lfc)  <- paste0(file_name, "_LFC")
  colnames(df_padj) <- paste0(file_name, "_padj")
  df_lfc$genes  <- rownames(data)
  df_padj$genes <- rownames(data)
  df_lfc  <- df_lfc  %>% filter(genes %in% jak)
  df_padj <- df_padj %>% filter(genes %in% jak)
  rownames(df_lfc)  <- df_lfc$genes;  df_lfc$genes  <- NULL
  rownames(df_padj) <- df_padj$genes; df_padj$genes <- NULL
  data_frames[[file_name]]  <- df_lfc
  padj_frames[[file_name]]  <- df_padj
}

all_row_names <- unique(unlist(lapply(data_frames, rownames)))
template_df   <- data.frame(row.names = all_row_names)
template_padj <- data.frame(row.names = all_row_names)

for (name in names(data_frames)) {
  df <- data_frames[[name]]
  template_df[[name]] <- df[match(all_row_names, rownames(df)), ]
  
  dfp <- padj_frames[[name]]
  template_padj[[name]] <- dfp[match(all_row_names, rownames(dfp)), ]
}

res  <- template_df
res_padj <- template_padj

#res_filtered <- res[rowSums(!is.na(res)) >= 5, ]
#res_padj_filtered <- res_padj[rownames(res_filtered), , drop=FALSE]

res_ordered <- res[, rownames(order_dep), drop = FALSE]
padj_ordered <- res_padj[, rownames(order_dep), drop = FALSE]

d <- as.matrix(res_ordered)
padj_mat <- as.matrix(padj_ordered)


stars <- ifelse(padj_mat < 0.05, "*", "")

minv <- -1
maxv <- 1
neutral_value <- 0
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue","lightblue"))(n = length(bk1)-1),
                "#e1e1e1","#e1e1e1",
                colorRampPalette(colors = c("tomato1","darkred"))(n = length(bk2)-1))

pheatmap(d, cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = bk, color = my_palette,
         annotation_col = order_dep,
         na_col = "#FFFFFF",
         display_numbers = stars, number_color = "black")

## segnare la signficatività a lato della heatmap 
## mut a sinistra e wt a destra

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut$genes <- rownames(mut)
mut <- mut %>% filter(genes %in% jak)
for (i in rownames(mut)) {
  if (mut[i, "padj"] < 0.05) {
    mut[i, "sign"] <- "sign"
  } else {
    mut[i, "sign"] <- "no"
  }
}
mut <- mut[, c("genes", "sign")]
mut$genes <- NULL
mut <- as.data.frame(t(mut))
mut$IFNL1 <- "NA"
mut <- as.data.frame(t(mut))
mut$genes <- rownames(mut)
mut <- mut[,c(2,1)]

wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE//WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt$genes <- rownames(wt)
wt <- wt %>% filter(genes %in% jak)
for (i in rownames(wt)) {
  if (wt[i, "padj"] < 0.05) {
    wt[i, "sign"] <- "sign"
  } else {
    wt[i, "sign"] <- "no"
  }
}

wt <- wt[, c("genes", "sign")]


col_fun <- colorRamp2(c(-1.5, 0, 1.5),
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
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80", "NA"="black"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)

draw(ht, merge_legend = TRUE)

analysis <- d
analysis <- as.data.frame(t(analysis))
analysis$sample <- rownames(analysis)
order_dep$sample <- rownames(order_dep)
analysis <- merge(analysis, order_dep, by="sample")
analysis$quartile <- NULL
wt <- analysis %>% filter(mut == "WT")
wt$mut <- NULL
rownames(wt) <- wt$sample
wt$sample <- NULL
wt <- as.data.frame(t(wt))
mut <- analysis %>% filter(mut == "MUT")
mut$mut <- NULL
rownames(mut) <- mut$sample
mut$sample <- NULL
mut <- as.data.frame(t(mut))

mut <- mut[rowSums(is.na(mut)) != ncol(mut), ]
#eliminate <- c("MIR221", "IL6", "IFNL1")
wt <- as.data.frame(t(wt))
wt$MIR221 <- NULL
wt$IL6 <- NULL
wt$IFNL1 <- NULL
wt <- as.data.frame(t(wt))
mut <- mut[rownames(wt), ]

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_jak <- data.frame(
  gene = names(results),
  value = unlist(results)
)

## adhesion
wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt <- wt %>% filter(p.adjust < 0.05)

ad <- wt[grep("adhesion", wt$Description),]
ad <- ad$geneID
ad <- as.character(ad)
ad <- unlist(strsplit(ad, "/"))
ad <- unique(ad)

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut <- mut %>% filter(p.adjust < 0.05)

ad2 <- mut[grep("adhesion", mut$Description),]
ad2 <- ad2$geneID
ad2 <- as.character(ad2)
ad2 <- unlist(strsplit(ad2, "/"))
ad2 <- unique(ad2)

ad <- unique(c(ad, ad2))

directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/logfc_singoli_N2.vs.NE_def/"

tsv_files <- list.files(directory_path, pattern = "\\.tsv$", full.names = TRUE)
data_frames <- list()

for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  data_frames[[file_name]] <- data
}

## tolgo CRC152 per FP

data_frames[["CRC0152"]] <- NULL

for (name in names(data_frames)) {
  df <- data_frames[[name]]
  df <- df[, 2, drop = FALSE]
  colnames(df) <- paste0(name, "_LFC")
  df$genes <- rownames(df)
  df <- df %>% filter(genes %in% ad)
  df$genes <- NULL
  data_frames[[name]] <- df
}

all_row_names <- unique(unlist(lapply(data_frames, rownames)))

template_df <- data.frame(row.names = all_row_names)

for (name in names(data_frames)) {
  df <- data_frames[[name]]
  template_df[[name]] <- df[match(all_row_names, rownames(df)), ]
}

res <- template_df

res_filtered <- res[rowSums(!is.na(res)) >= 15, ]
#res_filtered <- as.data.frame(t(res_filtered))

res_ordered <- res_filtered[,rownames(order_dep), drop = FALSE]

d <- res_ordered
pheatmap(d, cluster_rows = FALSE, cluster_cols = FALSE)
d <- as.matrix(d)

minv <- -5
maxv <- 5
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#e1e1e1", "#e1e1e1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette, annotation_col = order_dep, na_col = "#FFFFFF",
         fontsize_row = 5)


## extracellular matrix
wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt <- wt %>% filter(p.adjust < 0.05)

ecm <- wt[grep("extracellular matrix", wt$Description),]
ecm <- ecm$geneID
ecm <- as.character(ecm)
ecm <- unlist(strsplit(ecm, "/"))
ecm <- unique(ecm)

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut <- mut %>% filter(p.adjust < 0.05)

ecm2 <- mut[grep("extracellular matrix", mut$Description),]
ecm2 <- ecm2$geneID
ecm2 <- as.character(ecm2)
ecm2 <- unlist(strsplit(ecm2, "/"))
ecm2 <- unique(ecm2)

ecm <- unique(c(ecm, ecm2))

directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/logfc_singoli_N2.vs.NE_def/"

tsv_files <- list.files(directory_path, pattern = "\\.tsv$", full.names = TRUE)
data_frames <- list()
padj_frames <- list()

for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  df_lfc  <- data[, 2, drop = FALSE]
  df_padj <- data[, 6, drop = FALSE]
  colnames(df_lfc)  <- paste0(file_name, "_LFC")
  colnames(df_padj) <- paste0(file_name, "_padj")
  df_lfc$genes  <- rownames(data)
  df_padj$genes <- rownames(data)
  df_lfc  <- df_lfc  %>% filter(genes %in% ecm)
  df_padj <- df_padj %>% filter(genes %in% ecm)
  rownames(df_lfc)  <- df_lfc$genes;  df_lfc$genes  <- NULL
  rownames(df_padj) <- df_padj$genes; df_padj$genes <- NULL
  data_frames[[file_name]] <- df_lfc
  padj_frames[[file_name]] <- df_padj
}

all_row_names <- unique(unlist(lapply(data_frames, rownames)))
template_df   <- data.frame(row.names = all_row_names)
template_padj <- data.frame(row.names = all_row_names)

for (name in names(data_frames)) {
  df <- data_frames[[name]]
  template_df[[name]] <- df[match(all_row_names, rownames(df)), ]
  
  dfp <- padj_frames[[name]]
  template_padj[[name]] <- dfp[match(all_row_names, rownames(dfp)), ]
}

res <- template_df
res_padj <- template_padj

## togliere quei geni che non hanno lfc per almeno 15/30 campioni 
#res_filtered <- res[rowSums(!is.na(res)) >= 15, ] 
#res_filtered <- as.data.frame(t(res_filtered))
res_filtered <- res
res_padj_filtered <- res_padj[rownames(res_filtered), , drop=FALSE]

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

res_ordered  <- res_filtered[, rownames(order_dep), drop = FALSE]
padj_ordered <- res_padj_filtered[, rownames(order_dep), drop = FALSE]

d <- as.matrix(res_ordered)
padj_mat <- as.matrix(padj_ordered)

stars <- ifelse(padj_mat < 0.05, "*", "")

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
         annotation_col = order_dep, na_col = "#FFFFFF",
         display_numbers = stars, number_color = "black")


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
mut <- mut %>% filter(genes %in% ecm)
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
wt <- wt %>% filter(genes %in% ecm)
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

col_fun <- colorRamp2(c(-1.5, 0, 1.5),
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
  row_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(MUT = row_anno_mut,
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)

draw(ht, merge_legend = TRUE)

analysis <- d
analysis <- as.data.frame(t(analysis))
analysis$sample <- rownames(analysis)
order_dep$sample <- rownames(order_dep)
analysis <- merge(analysis, order_dep, by="sample")
analysis$quartile <- NULL
wt <- analysis %>% filter(mut == "WT")
wt$mut <- NULL
rownames(wt) <- wt$sample
wt$sample <- NULL
wt <- as.data.frame(t(wt))
mut <- analysis %>% filter(mut == "MUT")
mut$mut <- NULL
rownames(mut) <- mut$sample
mut$sample <- NULL
mut <- as.data.frame(t(mut))

#mut <- mut[rowSums(is.na(mut)) != ncol(mut), ]
#wt <- wt[!(rownames(wt) %in% "TMEM88B"), ]
mut <- mut[rownames(wt), ]

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_bcat_progress <- data.frame(
  gene = names(results),
  value = unlist(results)
)
