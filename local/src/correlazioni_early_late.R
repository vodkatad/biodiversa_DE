library(tidyverse)
library(ggrastr)

vsd_f <- snakemake@input[["vsd_file"]]
samples_f <- snakemake@input[["samples_file"]]
deg_f <- snakemake@input[["tsv_deg"]]
pdf1_p <- snakemake@output[["pdf1"]]
pdf2_p <- snakemake@output[["pdf2"]]
res <- snakemake@output[["tsv"]]

#vsd <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/early_late/vsd.tsv.gz", quote = "", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
vsd <- read.table(vsd_f, quote = "", sep = "\t",
                  header = TRUE, stringsAsFactors = FALSE)
#samples <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/early_late/samples_data", quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- read.table(samples_f, quote = "", sep = "\t",
                      header = TRUE, stringsAsFactors = FALSE)

# vsd <- as.data.frame(t(vsd))
# vsd$id <- rownames(vsd)
# merged <- merge(vsd, samples, by="id")
# new <- merged
# new$id <- NULL


### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
means <- apply(vsd, 1, mean)
med <- median(means)

de <- vsd[means > med,]

sds <- apply(de, 1, sd)

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

desd1 <- desd
#names(desd1) <- substr(names(desd1), 1, 12)
desd1 <- as.data.frame(t(desd1))
desd1$id <- rownames(desd1)
merged <- merge(desd1, samples, by="id")
new <- merged
new$model_passage <- paste0(new$model, "_", new$passage)
new$id <- NULL
#new$model <- NULL
new$passage <- NULL
rownames(new) <- new$model_passage
new$model_passage <- NULL


# creo una lista splittando per i modelli
list_df <- list()
list_df <- split(new, new$model)

# ## nuova colonna con modello e passaggio
# list_df <- lapply(list_df, function(new) {
#   new$model_passage <- rownames(new)
#   return(new)
# })

## sistemazione colonne
list_df <- lapply(list_df, function(new) {
  new$model <- NULL
  # new$passage <- NULL
  # rownames(new) <- new$model_passage
  # new$model_passage <- NULL
  new <- as.data.frame(t(new))
  return(new)
})

## cor test per la colonna 1 e 2 per tutti i df
cor_list <- list()
cor_list <- lapply(list_df, function(df) {
  cor <- cor.test(df[[1]], df[[2]])
  return(data.frame(cor = cor$estimate, pval = cor$p.value))
})

## metto il nome del df come colonna per non perdermi il sample
cor_list_with_names <- lapply(names(cor_list), function(df_name) {
  df <- cor_list[[df_name]]
  df$DataFrameName <- df_name
  return(df)
})

## creo un unico df
combined_df <- bind_rows(cor_list, .id = "DataFrameName")
rownames(combined_df) <- combined_df$DataFrameName
combined_df$DataFrameName <- NULL
combined_df <- combined_df[order(combined_df$cor),]

#deg_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/early_late/passage_cutoff0.05-early.vs.late.deseq2.tsv"
deg <- read.table(deg_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg <- deg %>% filter(padj < 0.05)
deg <- deg %>% filter(log2FoldChange > 0.5849625)
deg <- rownames(deg)

list_df <- lapply(list_df, function(new) {
  new$genes <- rownames(new)
  return(new)
})

specific_column <- "genes"

add_yes_no_column <- function(df, column_name, target_vector) {
  df %>%
    mutate(new_column = ifelse(.data[[column_name]] %in% target_vector, "yes", "no"))
}

list_df_deg <- lapply(list_df, add_yes_no_column, column_name = specific_column, target_vector = deg)

## funzione scatterplot
create_scatter_plot <- function(df) {
  ggplot(df, aes(x = df[[1]], y = df[[2]])) +
    rasterize(geom_point(size = 1, aes(color=df[[4]])), dpi=300) +
    geom_smooth(method = 'lm', size = 1) +
    xlab(names(df[1])) + ylab(names(df[2]))+
    guides(color = guide_legend(title = "DEG_gene"))+
    scale_color_manual(values = c("no" = "black", "yes" = "red"))
}
scatter_plot_list <- lapply(list_df_deg, create_scatter_plot)
pdf(file = pdf1_p)
print(scatter_plot_list[[1]])
dev.off()

pdf(file = pdf2_p)
print(scatter_plot_list[[16]])
dev.off()
write.table(combined_df, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)


