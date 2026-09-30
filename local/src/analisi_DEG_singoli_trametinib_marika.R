## check DEGS estremi trametinib marika

crc1139 <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/casi_singoli/results/CRC1139/TRAM_vs_NT/CRC1139_treat_cutoff0.05-TRAM.vs.NT.deseq2.tsv"
crc1139 <- read.table(crc1139, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

crc1139up <- crc1139 %>% filter(padj < 0.05)
crc1139up <- crc1139up %>% filter(log2FoldChange > 0.58)

crc1139down <- crc1139 %>% filter(padj < 0.05)
crc1139down <- crc1139down %>% filter(log2FoldChange < - 0.58)

crc1139$genes <- rownames(crc1139)

crc1620 <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/casi_singoli/results/CRC1620/TRAM_vs_NT/CRC1620_treat_cutoff0.05-TRAM.vs.NT.deseq2.tsv"
crc1620 <- read.table(crc1620, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

crc1620up <- crc1620 %>% filter(padj < 0.05)
crc1620up <- crc1620up %>% filter(log2FoldChange > 0.58)

crc1620down <- crc1620 %>% filter(padj < 0.05)
crc1620down <- crc1620down %>% filter(log2FoldChange < - 0.58)

crc1620$genes <- rownames(crc1620)

library(ggVennDiagram)

## UP

venn_list <- list(
  CRC1139 = rownames(crc1139up),
  CRC1620 = rownames(crc1620up)
)

p <- ggVennDiagram(venn_list)
p <- p + ggtitle("UP")
p +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(40, 60, 40, 60)  # aumenta soprattutto left/right
  )


## DOWN

venn_list <- list(
  CRC1139 = rownames(crc1139down),
  CRC1620 = rownames(crc1620down)
)

p <- ggVennDiagram(venn_list)
p <- p + ggtitle("DOWN")
p +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(40, 60, 40, 60)  # aumenta soprattutto left/right
  )

crc1139_df <- crc1139 %>%
  rename(
    log2FC_crc1139 = log2FoldChange,
    padj_crc1139 = padj
  )

crc1620_df <- crc1620 %>%
  rename(
    log2FC_crc1620 = log2FoldChange,
    padj_crc1620 = padj
  )

df <- merge(crc1139_df, crc1620_df, by = "genes")

df <- df %>%
  mutate(
    DEG_crc1139 = padj_crc1139 < 0.05 & abs(log2FC_crc1139) > 0.58,
    DEG_crc1620 = padj_crc1620 < 0.05 & abs(log2FC_crc1620) > 0.58
  ) %>%
  mutate(
    DEG_pattern = case_when(
      DEG_crc1139 & DEG_crc1620 ~ "both_models",
      DEG_crc1139 & !DEG_crc1620 ~ "crc1139_only",
      !DEG_crc1139 & DEG_crc1620 ~ "crc1620_only",
      TRUE ~ "none"
    )
  )

ggplot(df,
       aes(x = log2FC_crc1139,
           y = log2FC_crc1620,
           color = DEG_pattern)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    none = "grey70",
    crc1139_only = "#1f77b4",
    crc1620_only = "#d62728",
    both_models = "#2ca02c"
  )) +
  theme_minimal() + scale_y_continuous(limits = c(-8,8), expand = c(0,0)) + scale_x_continuous(limits = c(-8,8), expand = c(0,0))

ggplot(df, aes(log2FC_crc1139, log2FC_crc1620)) +
  rasterize(geom_point(alpha = 0.7), dpi=300) +
  theme_minimal() + scale_y_continuous(limits = c(-8,8), expand = c(0,0)) + scale_x_continuous(limits = c(-8,8), expand = c(0,0))

cor.test(df$log2FC_crc1139, df$log2FC_crc1620)

gsea_1139 <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/casi_singoli/results/CRC1139/TRAM_vs_NT/GSEA_results_CRC1139_H_treat_cutoff0.05-TRAM.vs.NT.tsv"
gsea_1139 <- read.table(gsea_1139, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_1139 <- gsea_1139[,c("ID", "NES", "p.adjust")]
colnames(gsea_1139) <- c("ID", "NES_1139", "p.adjust_1139")

gsea_1620 <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/casi_singoli/results/CRC1620/TRAM_vs_NT/GSEA_results_CRC1620_H_treat_cutoff0.05-TRAM.vs.NT.tsv"
gsea_1620 <- read.table(gsea_1620, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_1620 <- gsea_1620[,c("ID", "NES", "p.adjust")]
colnames(gsea_1620) <- c("ID", "NES_1620", "p.adjust_1620")

per_heatmap <- merge(gsea_1139, gsea_1620, by="ID")
rownames(per_heatmap) <- per_heatmap$ID
per_heatmap$ID <- NULL

nes_cols <- c("NES_1139", "NES_1620")
pval_cols <- c("p.adjust_1139","p.adjust_1620")

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
ord <- order(mat_nes[, "NES_1139"], decreasing = TRUE)


ht_nes <- Heatmap(mat_nes,
                  name = "NES",
                  col = colori,
                  #row_order=ord,
                  cluster_columns = FALSE,
                  cluster_rows = TRUE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
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
  grep("_1139$", colnames(mat_pval_cat), value = TRUE),
  grep("_1620$", colnames(mat_pval_cat), value = TRUE)
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


