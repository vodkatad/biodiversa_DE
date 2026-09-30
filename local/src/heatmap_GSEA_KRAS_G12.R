# mono <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/MRTX_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-MRTX.vs.NT.tsv"
# data_mono <- read.table(mono, header = TRUE, sep = "\t")

gsea_path<-"/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C.vs.NT.tsv"
data_gsea<-read.table(gsea_path,header = TRUE,sep='\t')

gsea_path_tripla<- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_T_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C_T.vs.NT.tsv"
gsea_path_tripla<-read.table(gsea_path_tripla,header = TRUE,sep='\t')

# gsea_tram <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/C_T_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-C_T.vs.NT.tsv"
# gsea_tram <- read.table(gsea_tram, header = TRUE, sep = "\t")

gsea_alp <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_B_vs_NT/GSEA_results_CRC0031_H_treat_cutoff0.05-M_C_B.vs.NT.tsv"
gsea_alp <- read.table(gsea_alp, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# m <- data_mono[,c('ID','NES','p.adjust')]
# colnames(m) <- paste0(colnames(m), "_MRTX")

mc<-data_gsea[,c('ID','NES','p.adjust')]
colnames(mc)<-paste0(colnames(mc),'_MRTX+CET')

# per_heatmap<-merge(m,mc,by="row.names")
# rownames(per_heatmap)<-per_heatmap$Row.names
# per_heatmap$Row.names<-NULL
# per_heatmap$ID_MRTX <- NULL
# per_heatmap$`ID_MRTX+CET` <- NULL

mct <- gsea_path_tripla[,c('ID','NES','p.adjust')]
colnames(mct) <- paste0(colnames(mct),'_MRTX+CET+TRAM')

per_heatmap<-merge(mc,mct,by="row.names")
rownames(per_heatmap)<-per_heatmap$Row.names
per_heatmap$Row.names<-NULL
per_heatmap$`ID_MRTX+CET+TRAM` <- NULL
per_heatmap$`ID_MRTX+CET` <- NULL

mca <- gsea_alp[,c('ID','NES','p.adjust')]
colnames(mca) <- paste0(colnames(mca),'_MRTX+CET+ALP')

per_heatmap <- merge(per_heatmap, mca, by="row.names")
rownames(per_heatmap) <- per_heatmap$Row.names
per_heatmap$Row.names <- NULL
per_heatmap$`ID_MRTX+CET+ALP` <- NULL

nes_cols <- c("NES_MRTX+CET", "NES_MRTX+CET+TRAM", "NES_MRTX+CET+ALP")
pval_cols <- c("p.adjust_MRTX+CET","p.adjust_MRTX+CET+TRAM", "p.adjust_MRTX+CET+ALP")

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
ord <- order(mat_nes[, "NES_MRTX+CET"], decreasing = TRUE)


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
                    labels = c(paste0("-", round(2, 2)), "0", paste0("+", round(2, 2))))
)
#print(ht_nes)

mat_pval<-as.matrix(mat_pval)
pval_th <- 0.05
mat_pval_cat <- ifelse(mat_pval < pval_th, "significativo", "non significativo")
mat_pval_cat[is.na(mat_pval)] <- "NA" 

mat_pval_cat <- as.matrix(mat_pval_cat)

ordered_columns <- c(
  grep("_MRTX\\+CET$", colnames(mat_pval_cat), value = TRUE),
  grep("_MRTX\\+CET\\+TRAM$", colnames(mat_pval_cat), value = TRUE),
  grep("_MRTX\\+CET\\+ALP$", colnames(mat_pval_cat), value = TRUE)
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
