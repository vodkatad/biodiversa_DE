library(tidyverse)

expr_input <- snakemake@input[["expr_input"]]
metadata_input <- snakemake@input[["metadata"]]
buoni <- snakemake@input[["goodOrBad"]]
expr_output <- snakemake@output[["expr_output"]]
metadata_output <- snakemake@output[["metadata_clu"]]
image <- snakemake@output[["script"]]

selection <- read.table(buoni, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote='')
#metadata_input <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/samples_datas"
metadata <- read.table(metadata_input, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote='')
#expr_input <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/vsd.tsv.gz"
#expr_df <- read.table(gzfile(expr_input), sep="\t", header=TRUE, stringsAsFactors=FALSE)
expr_df <- read.table(expr_input, sep="\t", header=TRUE, stringsAsFactors=FALSE)

good <- selection[selection$buoni, 'smodel']

lmo_good <- metadata[metadata$sample %in% good & metadata$type=="LMO_BASALE", ]

lmx_good <- metadata[metadata$sample %in% lmo_good$sample & metadata$type=="LMX_BASALE", ]
#lmh_good <- metadata[metadata$sample %in% lmo_good$sample & metadata$type=="LMH", ]

#has_h_has_x <- c(lmx_good$sample, lmh_good$sample)

#lmo_good_has_x <- lmo_good[lmo_good$sample %in% lmx_good$sample,]

wanted <- rbind(lmo_good, lmx_good)

expr_selected <- expr_df[,colnames(expr_df) %in% rownames(wanted)]

#class <- substr(colnames(expr_selected), 8, 10)
#table(class)
save.image(image)

### calcolo mediana della media e rimozione geni non espressi (expr media < mediana della media)
### procediamo solo con geni con sd alta
### NOTA: filtro fatto post selezione dei campioni quindi non coincidente al 100% con quello che avevamo
### per le correlazioni
# codice da vecchio filtro pre correlazioni che facevamo separatmente sulle matrici o/x
### da qui in poi , splittare e il tutto dovrebbe finire in una funzione per scegliere la sd
### delle tre classi
#######################################################################
means_expr <- apply(expr_selected, 1, mean)
med_expr <- median(means_expr)

### I obtain a data.table that contain only the genes  mean expression over the median
expr_filtered <- as.data.frame(expr_selected[means_expr > med_expr,])

### I want to calculate the standard deviation
sds_expr <- apply(expr_filtered, 1, sd) 

### now we keep the top 10% variable genes 
sds_expr <- sds_expr[order(-sds_expr)]
n_expr <- length(sds_expr)
keep_expr <- head(sds_expr, round(0.10*n_expr))
keep_genes_expr <- names(keep_expr)
expr_selected_fin <- expr_filtered[rownames(expr_filtered) %in% keep_genes_expr,]
###############################################################

write.table(expr_selected_fin, gzfile(expr_output), sep="\t", quote=FALSE)


write.table(wanted, metadata_output, sep="\t", quote=FALSE)