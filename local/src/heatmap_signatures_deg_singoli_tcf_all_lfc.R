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
bcat_go <- c("ZNRF3", "KREMEN1", "NKD1", "NOTUM", "HHEX", "PTK7", "LGR5", "LGR6", "APOE",  
             "NID1", "ADGRA2", "WNT8B", "DKK3", "DLX5", "HESX1", "DACT3", "FZD10", "TERT",  
             "BARX1", "FZD9", "RSPO4", "TMEM131L", "DISC1", "WNT10B", "FZD2", "GPRC5B", "LRP4",  
             "TPBGL", "DKKL1", "TMEM88B", "ARHGEF19")
bcat_progress <- c("LGR5", "EPHB3", "NKD1", "ZNRF3", "ASCL2", "GINS2",
                  "MYC", "SP5", "GINS3", "CFCA4", "TEAD4", "BMP4", "RNF43", 
                 "AXIN2", "CCND1", "DKK1")


res <- res %>% filter(genes %in% bcat_go)
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
orderdepcocomp <- order_dep
orderdepcocomp$quartile <- NULL
orderdepcocomp$mut <- NULL
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
mut <- mut %>% filter(genes %in% bcat_go)
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
wt <- wt %>% filter(genes %in% bcat_go)
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

# new order
cols_of_interest <- which(colnames(d) == "CRC1430") : which(colnames(d) == "CRC0161")
row_means_subset <- rowMeans(d[, cols_of_interest], na.rm = TRUE)
d <- d[order(row_means_subset, decreasing = FALSE), ]
row_anno_mut <- mut[match(rownames(d), mut$genes), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

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

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

mut_metagene <- as.data.frame(mut_metagene)
colnames(mut_metagene) <- c("metagene")
mut_metagene$type <- "MUT"
wt_metagene <- as.data.frame(wt_metagene)
colnames(wt_metagene) <- c("metagene")
wt_metagene$type <- "WT"

metagene_tot <- rbind(wt_metagene, mut_metagene)
metagene_tot$model <- rownames(metagene_tot)
metagene_tot <- merge(metagene_tot, orderdepcocomp, by="model")
metagene_plot <- metagene_tot
cor.test(metagene_tot$metagene, metagene_tot$co_comp_N2)

ggplot(metagene_tot, aes(x=metagene, y=co_comp_N2, color=type))+geom_point()+geom_smooth(method=lm) +geom_text_repel(aes(label = model), size = 3)

cor_results <- metagene_tot %>%
  group_by(type) %>%
  summarise(
    cor = cor(metagene, co_comp_N2, use = "complete.obs", method = "spearman"),
    p.value = cor.test(metagene, co_comp_N2)$p.value
  )

metagene_plot <- metagene_plot %>% filter(type == "WT")
metagene_plot$type <- NULL
ggplot(metagene_plot, aes(x=metagene, y=co_comp_N2))+geom_point()+geom_smooth(method=lm) +geom_text_repel(aes(label = model), size = 3)
cor.test(metagene_plot$metagene, metagene_plot$co_comp_N2)
results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_bcat <- data.frame(
  gene = names(results),
  value = unlist(results)
)

results_bcat_sing <- results_bcat %>% filter(value<0.05)

## parte heatmap jak

res <- lfc_res

jak <- c("PTK2B", "IL23R", "IL6R", "CCL5", "IL7R", "TNF", "MIR221", "GHR", "IFNL1", "IL6") 

res <- res %>% filter(genes %in% jak)
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
mut <- mut %>% filter(genes %in% jak)
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
wt <- wt %>% filter(genes %in% jak)
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

# new order
wt_cols <- colnames(d)[13:30]

co_comp_WT <- orderdepcocomp[wt_cols, "co_comp_N2"]

gene_cor <- apply(d, 1, function(x) {
  vals <- x[wt_cols]
  if (all(is.na(vals))) return(NA)
  cor(vals, co_comp_WT)
})


d <- d[order(gene_cor, decreasing = FALSE), ]

cor_df <- data.frame(
  gene = rownames(d),
  cor_WT = gene_cor[match(rownames(d), names(gene_cor))]
)

row_anno_mut <- mut[match(rownames(d), mut$genes), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

top_anno <- HeatmapAnnotation(df = order_dep[,c(1,2)])

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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_jak <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
tnfmut <- mut %>% filter(rownames(mut)=="TNF")
tnfmut <- as.numeric(tnfmut[1,])
tnfwt <- wt %>% filter(rownames(wt)=="TNF")
tnfwt <- as.numeric(tnfwt[1,])
wilcox.test(tnfmut, tnfwt)
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)

wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

mut_metagene <- as.data.frame(mut_metagene)
colnames(mut_metagene) <- c("metagene")
mut_metagene$type <- "MUT"
wt_metagene <- as.data.frame(wt_metagene)
colnames(wt_metagene) <- c("metagene")
wt_metagene$type <- "WT"

metagene_tot <- rbind(wt_metagene, mut_metagene)
metagene_tot$model <- rownames(metagene_tot)
metagene_tot <- merge(metagene_tot, orderdepcocomp, by="model")
cor.test(metagene_tot$metagene, metagene_tot$co_comp_N2)

ggplot(metagene_tot, aes(x=metagene, y=co_comp_N2, color=type))+geom_point()+geom_smooth(method=lm) +geom_text_repel(aes(label = model), size = 3)

cor_results <- metagene_tot %>%
  group_by(type) %>%
  summarise(
    cor = cor(metagene, co_comp_N2, use = "complete.obs", method = "spearman"),
    p.value = cor.test(metagene, co_comp_N2)$p.value
  )

## citokyn mediated pathway
res <- lfc_res

#jak <- c("PTK2B", "IL23R", "IL6R", "CCL5", "IL7R", "TNF", "MIR221", "GHR", "IFNL1", "IL6") 
cit <- c(
 "PADI2",    "IL4R"  ,   "DUOX2",    "CEACAM1",  "CD24"  ,   "IL22RA1",  "MAPK3" ,   "PTK2B" ,   "OAS1",    
 "APPL2" ,   "CCRL2" ,   "ACKR2" ,   "IFI27",    "F2RL1",    "OASL",     "STAT2" ,   "IRF7",     "TMSB4X",
 "IL2RG" ,   "IFNLR1",   "SLC1A1" ,  "NFKBIZ" ,  "CCL20"  ,  "CLCF1" ,   "GPR17" ,   "RPS6KA5",  "CXCL11" , 
 "EDN1" ,    "IL1R2" ,   "MX1" ,     "HPX" ,     "ISG15" ,   "TFF2",     "OAS2"  ,   "BIRC3" ,   "IL23R" ,  
 "IL33" ,    "APOA1" ,   "IL1RN",    "IL6R" ,    "CCL5" ,    "TNFRSF14", "MMP12" ,   "SPHK1" ,   "CCL14" ,  
 "CCL22",    "ZBP1" ,    "IL36RN" ,  "IL1B" ,    "IL18R1" ,  "NR1H4" ,   "CARD16" ,  "USP18" ,   "CXCL1" ,  
 "CASP1" ,   "IL7R" ,    "PPBP",     "CSF1" ,    "CXCL10" ,  "CLDN18" ,  "MIR21" ,   "DUOX1" ,   "TNF"  ,   
 "CXCL8" ,   "IL17REL" , "EDN2" ,    "CCR9" ,    "GHR",      "NLRP6" ,   "IL13RA2" , "INHBA" ,   "IFNL1" ,  
 "XCR1" ,    "IL1RL1" ,  "LIFR" ,    "CCL16",    "CSF2RB" ,  "IL1A"  ,   "IL6"  ,    "CCR8"  ,   "EBI3"  ,  
 "TICAM2",   "IL2RB" ,   "IL36B", "CNTFR", "CXCL6", "ILF1F10", "TNFSF11") 

res <- res %>% filter(genes %in% cit)
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

minv <- -2
maxv <- 2
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
mut <- mut %>% filter(genes %in% cit)
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
wt <- wt %>% filter(genes %in% cit)
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

col_fun <- colorRamp2(c(-2, 0, 2),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_cit <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

### heatmap ecm

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

res <- lfc_res

res <- res %>% filter(genes %in% ecm)
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

minv <- -3
maxv <- 3
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_ecm <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

## gruppo MHC signature solo wt
mhc <- c("B2M", "HLA-E", "HLA-F", "HLA-C", "HLA-B", "HLA-G", "HLA-DRB1", "ULBP2", "SERPINE2",
        "SLC6A4", "TAC1", "P2RX1", "SEMG1", "CHP1", "CAPN2", "FABP3", "NR1H4", "ADGRF5", 
        "IL1R2", "IL1RN", "ZBP1", "MIR21", "IL6", "S100A9", "ALOX5", "MDK", "TNF", "FFAR2", "S100A8") 

res <- lfc_res

res <- res %>% filter(genes %in% mhc)
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

minv <- -3
maxv <- 3
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
mut <- mut %>% filter(genes %in% mhc)
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
wt <- wt %>% filter(genes %in% mhc)
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_ecm <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")


## goblet da liste geni

res <- lfc_res

gob <- c("SPDEF", "TFF3", "CLCA1", "IL33", "MUC2", "MUC4",
         "FCGBP", "ZG16", "LGALS2", "AGR2", "SYTL2", "FER1L6",
         "MALAT1", "MUC13", "CREB3L1", "PGHR1", "GUCA2A", "KRT20")

gob <- c("SPDEF", "ATOH1", "CLCA1", "TFF3", "BEST2", "SPINK4", "REP15", "CREB3L1",
         "MUC2", "FCGBP", "ZG16", "BCAS1", "AGR2", "SYTL2", "MUC13",
         "KLF4", "REG4", "B4GALNT2", "MALAT1", "LGALS2", "FER1L6")

res <- res %>% filter(genes %in% gob)
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

minv <- -2
maxv <- 2
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
mut <- mut %>% filter(genes %in% gob)
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
wt <- wt %>% filter(genes %in% gob)
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_gob <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")


## geni ricavati dalle private di tcf
res <- lfc_res

wt <- "/home/egrassi/onlywt_genes_go_00005.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gowt <- wt$gene

res <- res %>% filter(genes %in% gowt)
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

minv <- -3
maxv <- 3
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
mut <- mut %>% filter(genes %in% gowt)
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
wt <- wt %>% filter(genes %in% gowt)
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_gowt <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

wt_cols <- colnames(d)[13:30]

co_comp_WT <- orderdepcocomp[wt_cols, "co_comp_N2"]

gene_cor <- apply(d, 1, function(x) {
  vals <- x[wt_cols]
  if (all(is.na(vals))) return(NA)
  cor(vals, co_comp_WT)
})


d <- d[order(gene_cor, decreasing = FALSE), ]

cor_df <- data.frame(
  gene = rownames(d),
  cor_WT = gene_cor[match(rownames(d), names(gene_cor))]
)

row_anno_mut <- mut[match(rownames(d), rownames(mut)), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

top_anno <- HeatmapAnnotation(df = order_dep[,c(1,2)])

ht <- Heatmap(
  as.matrix(d),
  name = "LFC",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  na_col = "#FFFFFF",
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 3),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(MUT = row_anno_mut,
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)
draw(ht, merge_legend = TRUE)
write.xlsx(cor_df, file="correlazione_dipendenza_geni_go_wt_private.xlsx")

## ligandi egf
res <- lfc_res

egf <- c("AREG", "EREG", "EGF", "HBEGF", "BTC", "TGFA")

res <- res %>% filter(genes %in% egf)
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

minv <- -3
maxv <- 3
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
mut <- mut %>% filter(genes %in% egf)
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
wt <- wt %>% filter(genes %in% egf)
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_egf <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

wt_cols <- colnames(d)[13:30]

co_comp_WT <- orderdepcocomp[wt_cols, "co_comp_N2"]

gene_cor <- apply(d, 1, function(x) {
  vals <- x[wt_cols]
  if (all(is.na(vals))) return(NA)
  cor(vals, co_comp_WT)
})


d <- d[order(gene_cor, decreasing = FALSE), ]

cor_df <- data.frame(
  gene = rownames(d),
  cor_WT = gene_cor[match(rownames(d), names(gene_cor))]
)

row_anno_mut <- mut[match(rownames(d), mut$genes), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

top_anno <- HeatmapAnnotation(df = order_dep[,c(1,2)])

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

### selezionate morte cellulare

res <- lfc_res

go_private <- "/home/egrassi/onlywt.tsv"
go <- read.table(go_private, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
morte <- c("response to oxidative stress", "response to hydroperoxide", "pyroptosis", "positive regulation of programmed cell death")
go <- go %>% filter(desc %in% morte)          
morte <- go$genes
morte <- as.character(morte)
morte <- unlist(strsplit(morte, "/"))
morte <- unique(morte)

res <- res %>% filter(genes %in% morte)
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

minv <- -3
maxv <- 3
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
mut <- mut %>% filter(genes %in% morte)
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
wt <- wt %>% filter(genes %in% morte)
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

col_fun <- colorRamp2(c(-3, 0, 3),
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

results <- list()

for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), na.action=na.exclude)$p.value
}

results_morte <- data.frame(
  gene = names(results),
  value = unlist(results)
)

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
wilcox.test(mut_metagene, wt_metagene)

boxplot(list(MUT = mut_metagene, WT = wt_metagene),
        ylab = "mean_logFC")

wt_cols <- colnames(d)[13:30]

co_comp_WT <- orderdepcocomp[wt_cols, "co_comp_N2"]

gene_cor <- apply(d, 1, function(x) {
  vals <- x[wt_cols]
  if (all(is.na(vals))) return(NA)
  cor(vals, co_comp_WT)
})


d <- d[order(gene_cor, decreasing = FALSE), ]

cor_df_morte <- data.frame(
  gene = rownames(d),
  cor_WT = gene_cor[match(rownames(d), names(gene_cor))]
)

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut$genes <- rownames(mut)
mut <- mut %>% filter(genes %in% morte)
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
wt <- wt %>% filter(genes %in% morte)
for (i in rownames(wt)) {
  if (wt[i, "padj"] < 0.05) {
    wt[i, "sign"] <- "sign"
  } else {
    wt[i, "sign"] <- "no"
  }
}

wt <- wt[, c("genes", "sign")]
row_anno_mut <- mut[match(rownames(d), mut$genes), "sign"]
row_anno_wt  <- wt[match(rownames(d), wt$genes), "sign"]

top_anno <- HeatmapAnnotation(df = order_dep[,c(1,2)])

ht <- Heatmap(
  as.matrix(d),
  name = "LFC",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  na_col = "#FFFFFF",
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 3),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(MUT = row_anno_mut,
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)
draw(ht, merge_legend = TRUE)
write.xlsx(cor_df, file="correlazione_dipendenza_geni_go_wt_private.xlsx")