library(tidyverse)
library(WriteXLS)

# we have 6 GO and 3 GSEA
list_enrich_results <- snakemake@input
go_gsea <- snakemake@output[["sign_go_gsea_xls"]] 
tsv <- snakemake@output[["merged_tsv"]] 
vs <- c('LMO-LMH','LMX-LMH','LMO-LMX')
vs <- c(vs, vs, vs)


load_preprocess_go <- function(filename, vs_name) {
  df <- read.table(filename, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  #df <- subset(df, p.adjust < 0.05) # rimuoverei il filtro sul pvalue per vedere gli interi risultati, in effetti in GSEA non lo fai e mi sembra in linea
  df$conc <- paste(df$ID, df$Description)
  col_order <- c("conc", 'pvalue')
  df <- df[, col_order]
  df <- df %>% remove_rownames %>% column_to_rownames(var="conc")
  names(df)[names(df) == 'pvalue'] <- paste0("pvalue", "_", vs_name)
  # We plug in positive or negative log(padj) depending on the direction 
  # of the genes
  # padj 0.002 -> log10 -> -3
  # so for up we do -(-3) = 3
  # for down the opposite
  return(df)
}

load_preprocess_gsea <- function(filename, vs_name) {
  df <- read.table(filename, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df$conc <- paste(df$ID, df$Description)
  # df$ID <- NULL
  # df$Description <- NULL
  # df$setSize <- NULL
  # df$enrichmentScore <- NULL
  # df$pvalue <- NULL
  # df$rank <- NULL
  # df$leading_edge <- NULL
  # df$core_enrichment <- NULL
  col_order <- c("conc", 'pvalue', "NES")
  df <- df[, col_order]
  df <- df %>% remove_rownames %>% column_to_rownames(var="conc")
  names(df)[names(df) == 'pvalue'] <- paste0("pvalue", "_", vs_name)
  names(df)[names(df) == 'NES'] <- paste0("NES", "_", vs_name)
  #colnames(df) <- c(vs_name, 'NES')
  #df <- df[df[,vs_name] < 0.05, , drop = FALSE]
  return(df)
}


for (i in seq(1,3)) {
  df_preproc <- load_preprocess_go(list_enrich_results[[i]], vs[i])
  if (i == 1) {
    go_merged_up <- df_preproc
  } else {
    go_merged_up <- merge(go_merged_up, df_preproc, by="row.names", all=TRUE)
    go_merged_up <- go_merged_up %>% remove_rownames %>% column_to_rownames(var="Row.names")
  }
}

for (i in seq(4,6)) {
  df_preproc <- load_preprocess_go(list_enrich_results[[i]], vs[i])
  if (i == 4) {
    go_merged_down <- df_preproc
  } else {
    go_merged_down <- merge(go_merged_down, df_preproc, by="row.names", all=TRUE) 
    go_merged_down <- go_merged_down %>% remove_rownames %>% column_to_rownames(var="Row.names")
  }
}

go_merged_up$sign <- "up"
go_merged_down$sign <- "down"
go_rbinded <- rbind(go_merged_up, go_merged_down)
#go_rbinded[is.na(go_rbinded)] <- 0.00

#write.table(go_rbinded, file = go, col.names = TRUE, row.names = TRUE,
            #sep = "\t", quote = FALSE)

# for (i in seq(7,9)) { ...same for gsea with ad hoc function ... }
for (i in seq(7,9)) {
  df_preproc <- load_preprocess_gsea(list_enrich_results[[i]], vs[i])
  if (i == 7) {
    gsea_merged <- df_preproc
  } else {
    gsea_merged <- merge(gsea_merged, df_preproc, by="row.names", all = TRUE)
    gsea_merged <- gsea_merged %>% remove_rownames %>% column_to_rownames(var="Row.names")
  }
}
go_rbinded <- tibble::rownames_to_column(go_rbinded, "Description")
gsea_merged <- tibble::rownames_to_column(gsea_merged, "Description")
#gsea_merged[is.na(gsea_merged)] <- 0.00

gsea_merged$typenrich <- 'GSEA'
go_rbinded$typenrich <- 'GO'


# We need to adjust all together pvalues for the three comparisons (X-O-H) and GSEA/GO
# so we put them together in a single vector, correct it, and then reassing the adjusted
# pvalues back to three columns with the same order as the pvalue ones.
merged <- dplyr::bind_rows(go_rbinded, gsea_merged)

col_pvals <- grep('pvalue_', colnames(merged))
for (i in seq(1, length(col_pvals))) {
  if (i == 1 ){
    all_pvals <- merged[,col_pvals[i]]    
  } else {
    all_pvals <- c(all_pvals, merged[,col_pvals[i]])
  }
}

all_padj <- p.adjust(all_pvals, method="BH")

df_padj <- as.data.frame(matrix(data=all_padj, nrow=nrow(merged), ncol=length(col_pvals), byrow=FALSE))
#df_pval <- as.data.frame(matrix(data=all_pvals, nrow=nrow(merged), ncol=length(col_pvals), byrow=FALSE))

col_padj <- gsub('pvalue_', 'padj_', colnames(merged)[col_pvals], fixed=TRUE)
colnames(df_padj) <- col_padj
merged <- cbind(merged, df_padj)

merged <- merged[order(merged$`pvalue_LMO-LMX`),]

#merged[is.na(merged)] <- "NA"  # non li abbiamo gia` convertiti prima a 0? Teoricamente qui se teniamo sempre tutti i termini non dovremmo averne di mancanti, non avendo nessun filtro sul pvalue.
write.table(merged, file = tsv, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

openxlsx::write.xlsx(merged, file = go_gsea, sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

