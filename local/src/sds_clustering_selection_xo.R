#setwd("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/")
library(tidyverse)

snakemakeorig <- snakemake
input_data <- snakemake@input[["script"]]

load(input_data)

expr_output <- snakemakeorig@output[["expr_output"]]
metadata_output <- snakemakeorig@output[["metadata_clu"]]

### calcolo mediana della media e rimozione geni non espressi (expr media < mediana della media)
### procediamo solo con geni con sd alta
### NOTA: filtro fatto post selezione dei campioni quindi non coincidente al 100% con quello che avevamo
### per le correlazioni
# codice da vecchio filtro pre correlazioni che facevamo separatmente sulle matrici o/x
### da qui in poi , splittare e il tutto dovrebbe finire in una funzione per scegliere la sd
### delle tre classi
### prendere lmx good, lmh good e lmohasxh
#######################################################################


expr_lmo <- expr_selected
expr_lmx <- expr_selected

expr_lmx <- expr_lmx[, colnames(expr_lmx) %in% rownames(lmx_good)]
expr_lmo <- expr_lmo[, colnames(expr_lmo) %in% rownames(lmo_good_has_xh)]

sd_analysis <- function(dataset) {
  means_expr <- apply(dataset, 1, mean)
  med_expr <- median(means_expr)
  ### I obtain a data.table that contain only the genes  mean expression over the median
  expr_filtered <- as.data.frame(dataset[means_expr > med_expr,])
  ### I want to calculate the standard deviation
  sds_expr <- apply(expr_filtered, 1, sd) 
  ### now we keep the top 10% variable genes 
  sds_expr <- sds_expr[order(-sds_expr)]
  n_expr <- length(sds_expr)
  keep_expr <- head(sds_expr, round(0.10*n_expr))
  keep_genes_expr <- names(keep_expr)
  expr_selected_fin <- expr_filtered[rownames(expr_filtered) %in% keep_genes_expr,]
  return(expr_selected_fin)
}

sds_lmx <- sd_analysis(expr_lmx)
sds_lmx <- tibble::rownames_to_column(sds_lmx, "genes")
sds_lmo <- sd_analysis(expr_lmo)
sds_lmo <- tibble::rownames_to_column(sds_lmo, "genes")  

sds_merged <- merge(sds_lmo, sds_lmx, by = "genes")

sds_merged <- sds_merged %>% remove_rownames %>% column_to_rownames(var = "genes")

###############################################################

write.table(sds_merged, gzfile(expr_output), sep="\t", quote=FALSE)


write.table(wanted, metadata_output, sep="\t", quote=FALSE)