## mettere fpkm
library(openxlsx)

geni <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/local/share/data/TCF7L2/geni_expr.xlsx"
geni <- read_xlsx(geni, col_names = FALSE)
colnames(geni) <- c("genes")
genes <- geni$genes

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2/fpkm.tsv.gz"
vsd <- read.table(vsd)
vsd$geni <- rownames(vsd)
vsd <- vsd %>% filter(geni %in% genes)
vsd$geni <- NULL
vsd <- as.data.frame(t(vsd))
vsd$model <- rownames(vsd)
vsd$geno <- substr(vsd$model, 9, 11)
vsd <- vsd %>% filter(geno %in% c("NKO", "CAS"))
vsd$model <- NULL
vsd$geno <- NULL
vsd <- log(vsd+1)
#vsd <- as.data.frame(t(vsd))

pheatmap(vsd, cluster_rows = FALSE, cluster_cols = FALSE)

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/general/fpkm.tsv.gz"
vsd <- read.table(vsd)
vsd$geni <- rownames(vsd)
vsd <- vsd %>% filter(geni %in% genes)
#> setdiff(genes, rownames(vsd))
#[1] "ASLC2"
vsd$geni <- NULL
vsd <- as.data.frame(t(vsd))
vsd$model <- substr(rownames(vsd), 1, 10)
vsd$replica <- substr(rownames(vsd), 12, 13)
vsd$model <- gsub("_", ".", vsd$model)
vsd$all <- paste0(vsd$model, "_", vsd$replica)
rownames(vsd) <- vsd$all
vsd$model <- NULL
vsd$replica <- NULL
vsd$all <- NULL
vsd <- log(vsd+1)
#vsd <- as.data.frame(t(vsd))
vsd$cases <- rownames(vsd)
vsd <- vsd[, c(ncol(vsd), 1:(ncol(vsd)-1))]

write.xlsx(vsd, file="/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/pheatmap_geneal.xlsx")

pheatmap(vsd, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 5)
