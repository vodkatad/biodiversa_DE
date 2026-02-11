library(dplyr)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(circlize)

# cartella con tutti gli excel
# /home/mferri/excel_geni_tcf_bcat

# fpkm singoli
directory_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/fpkm_singoli_N2.vs.NE_def"
# co_comp
order_dep_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx"
# DEGs globali MUT e WT
mut_path <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
wt_path  <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE//WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"

# lista geni
#gene_list <- c("PTK2B", "IL23R", "IL6R", "CCL5", "IL7R", "TNF", "MIR221", "GHR", "IFNL1", "IL6")
# go_private <- "/home/egrassi/onlywt.tsv"
# go <- read.table(go_private, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# morte <- c("pyroptosis", "positive regulation of programmed cell death")
# go <- go %>% filter(desc %in% morte)
# morte <- go$genes
# morte <- as.character(morte)
# morte <- unlist(strsplit(morte, "/"))
# morte <- unique(morte)
# gene_list <- morte
# wt <- "/home/egrassi/onlywt_genes_go_00005.tsv"
# wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# gene_list <- wt$gene
go_private <- "/home/egrassi/onlywt.tsv"
go <- read.table(go_private, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
morte <- c("response to oxidative stress", "response to hydroperoxide", "regulation of nitric oxide metabolic process",
           "regulation of nitric oxide biosynthetic process", "nitric oxide metabolic process", "nitric oxide biosynthetic process",
           "pyroptosis", "positive regulation of programmed cell death")
go <- go %>% filter(desc %in% morte)
morte <- go$genes
morte <- as.character(morte)
morte <- unlist(strsplit(morte, "/"))
morte <- unique(morte)
gene_list <- morte
# gene_list <- c("SPHK1","SPRR2A","NOS2","IL6","HLA-E","NLRP6",
#                "DHX58","PGLYRP4","HSPB1","CCL16","CSF2RB","LGALS9",
#                "S100A9","DDX60","IL7R","TRAV27","IGHA1","OASL",
#                "NODAL","IRF7","IL2RG","THBS1","IL4R","NEAT1")
# gene_list <- c("TMIGD1","NOS2","SPHK1","IL6","RCAN1")
# gene_list <- c("NOS2", "IL6","FANK1","NLRP6","LGALS9","S100A9",
#                "GADD45B","IFI27","SMPD1","BMF","THBS1","PHLDA3","CDKN1A")
gene_list <- c("SPHK1","SPRR2A","NOS2","IL6","HLA-E",
               "NLRP6","DHX58","PGLYRP4","HSPB1","CCL16",
               "CSF2RB","LGALS9","S100A9","DDX60","IL7R","TRAV27","IGHA1",
               "OASL","NODAL","IRF7","IL2RG","THBS1","IL4R","NEAT1","RCAN1","TMIGD1")

## calcolo lfc per tutti i singoli
tsv_files <- list.files(directory_path, pattern = "\\.gz$", full.names = TRUE)
data_list <- lapply(tsv_files, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  return(df)
})
names(data_list) <- gsub("_fpkm\\.tsv\\.gz$", "", basename(tsv_files))

colnames(data_list[["CRC0542"]]) <- c("CRC0542_NE_R1", "CRC0542_N2_R1", "CRC0542_NE_R2", "CRC0542_N2_R2")

calc_logFC_df <- function(df, sample_name, PC = 1) {
  cols_N2 <- grep("N2", colnames(df), value = TRUE)
  cols_NE <- grep("NE", colnames(df), value = TRUE)
  logFC_matrix <- mapply(function(n2, ne) log2((df[[n2]] + PC) / (df[[ne]] + PC)), cols_N2, cols_NE)
  logFC <- apply(logFC_matrix, 1, median)
  data.frame(row.names = rownames(df), logFC)
}

logFC_list <- mapply(calc_logFC_df, data_list, names(data_list), SIMPLIFY = FALSE)
names(logFC_list) <- gsub("_fpkm\\.tsv$", "", names(logFC_list))

all_row_names <- unique(unlist(lapply(logFC_list, rownames)))
lfc_res <- data.frame(row.names = all_row_names)
for (name in names(logFC_list)) {
  lfc_res[[paste0(name, "_LFC")]] <- logFC_list[[name]][match(all_row_names, rownames(logFC_list[[name]])), 1]
}
lfc_res$genes <- rownames(lfc_res)

res <- lfc_res %>% filter(genes %in% gene_list)
res$genes <- NULL
colnames(res) <- gsub("_LFC", "", colnames(res))

## annotation col 
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
order_dep <- order_dep %>% filter(!model == c("CRC0743", "CRC0152")) 
orderdepcocomp <- order_dep
orderdepcocomp$quartile <- NULL
order_dep$co_comp_N2 <- NULL 
order_dep$model <- NULL 
order_dep <- order_dep[order(order_dep$mut),] 
#order_dep$quartile <- NULL 
order_depnoq <- order_dep 
order_depnoq$quartile <- NULL

res_ordered  <- res[, rownames(order_dep), drop = FALSE]

d <- as.matrix(res_ordered)

## cor wt dipendenza
wt_cols <- colnames(d)[14:30]

co_comp_WT <- orderdepcocomp[wt_cols, "co_comp_N2"]
gene_cor <- apply(d, 1, function(x) {
  vals <- x[wt_cols]
  if (all(is.na(vals))) return(NA)
  cor(vals, co_comp_WT)
})
d <- d[order(gene_cor, decreasing = FALSE), ]

#annotazione wt e mut significativi
mut_anno <- read.table(mut_path, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut_anno$genes <- rownames(mut_anno)
mut_anno <- mut_anno %>% filter(genes %in% gene_list)
for (i in rownames(mut_anno)) {
  if (mut_anno[i, "padj"] < 0.05) {
    mut_anno[i, "sign"] <- "sign"
  } else {
    mut_anno[i, "sign"] <- "no"
  }
}
mut_anno <- mut_anno[, c("genes", "sign")]

wt_anno <- read.table(wt_path, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt_anno$genes <- rownames(wt_anno)
wt_anno <- wt_anno %>% filter(genes %in% gene_list)
for (i in rownames(wt_anno)) {
  if (wt_anno[i, "padj"] < 0.05) {
    wt_anno[i, "sign"] <- "sign"
  } else {
    wt_anno[i, "sign"] <- "no"
  }
}

wt_anno <- wt_anno[, c("genes", "sign")]

## heatmap
col_fun <- colorRamp2(c(-3, 0, 3),
                      c("darkblue", "#e1e1e1", "darkred"))

row_anno_mut <- mut_anno[match(rownames(d), mut_anno$genes), "sign"]
row_anno_wt  <- wt_anno[match(rownames(d), wt_anno$genes), "sign"]

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
  row_names_gp = gpar(fontsize = 6),
  top_annotation = top_anno,
  left_annotation = rowAnnotation(MUT = row_anno_mut,
                                  col = list(MUT = c("sign" = "tomato", "no" = "grey80"))),
  right_annotation = rowAnnotation(WT = row_anno_wt,
                                   col = list(WT = c("sign" = "steelblue", "no" = "grey80")))
)
draw(ht)

## boxplot
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

wt2 <- wt
mut2 <- mut
mut_metagene <- colMeans(mut2, na.rm = TRUE)
wt_metagene  <- colMeans(wt2,  na.rm = TRUE)
w <- wilcox.test(mut_metagene, wt_metagene)$p.value

load("vettore_metagene_bcat.rda")
boxplot(list(BCAT = mutagene_tot, WT = wt_metagene, MUT = mut_metagene),
        ylab = "Median_logFC", ylim =c(-0.2,1.5))
w2 <- wilcox.test(mutagene_tot,wt_metagene)$p.value
w3 <- wilcox.test(mutagene_tot, mut_metagene)$p.value

# correlazioni gene–dipendenza
cor_df <- data.frame(gene = names(gene_cor), cor_WT = gene_cor)
cor_df <- cor_df[order(cor_df$cor_WT),]
write.xlsx(cor_df, "morteestressetuttoquanto_correlation_dipendenza.xlsx", rowNames = FALSE)

# Wilcoxon test gene per gene
results <- list()
for (i in rownames(mut)) {
  results[[i]] <- wilcox.test(unlist(mut[i,]), unlist(wt[i,]), side="less", na.action=na.exclude)$p.value
}
results <- data.frame(
  gene = names(results),
  value = unlist(results)
)
results <- results[order(results$value),]
write.xlsx(results, "wilcoxon_morteestressetuttoquanto_pvalues.xlsx", rowNames = FALSE)

scatter <- as.data.frame(t(d))
scatter$row_mean <- rowMeans(scatter, na.rm = TRUE)
scatter$model <- rownames(scatter)
scatter <- scatter[,c("model", "row_mean")]
scatter <- merge(scatter, orderdepcocomp, by="model")

ggplot(scatter, aes(x=row_mean, y=co_comp_N2, color=mut))+geom_point()+geom_smooth(method = lm)+geom_text_repel(aes(label = model))+xlab("Mean_LFC_DEG_model")
