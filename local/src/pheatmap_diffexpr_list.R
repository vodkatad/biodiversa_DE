### heatmap shutte

library(pheatmap)
library(tidyverse)

# we have 6 GO and 3 GSEA
list_enrich_results <- snakemake@input
plot <- snakemake@output[["hm"]]
go_gsea <- snakemake@output[["sign_go_gsea"]]
go <- snakemake@output[["sign_go"]]
plot_20 <- snakemake@output[["hm_20"]]
vs <- c('LMO-LMH','LMX-LMH','LMO-LMX')
vs <- c(vs, vs, vs)

save.image('pippo.Rdata')

load_preprocess_go <- function(filename, vs_name) {
	df <- read.table(filename, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	df$X <- paste(df$ID, df$Description)
	df$ID <- NULL
	df$Description <- NULL
	df$GeneRatio <- NULL
	df$BgRatio <- NULL
	df$pvalue <- NULL
	df$qvalue <- NULL
	df$geneID <- NULL
	df$Count <- NULL
	col_order <- c("X", "p.adjust")
	df <- df[, col_order]
	df <- df %>% remove_rownames %>% column_to_rownames(var="X")
	colnames(df) <- c(vs_name)
	df <- df[df[,vs_name] < 0.05, , drop = FALSE]
	# We plug in positive or negative log(padj) depending on the direction 
	# of the genes
	# padj 0.002 -> log10 -> -3
	# so for up we do -(-3) = 3
	# for down the opposite
	if (grepl('up', filename, fixed=TRUE)) {
		df[,vs_name] <- -log10(df[, vs_name]) 
	} else {
		df[,vs_name] <- log10(df[, vs_name])
	}
	return(df)
}

load_preprocess_gsea <- function(filename, vs_name) {
	df <- read.table(filename, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	df$X <- paste(df$ID, df$Description)
	# df$ID <- NULL
	# df$Description <- NULL
	# df$setSize <- NULL
	# df$enrichmentScore <- NULL
	# df$pvalue <- NULL
	# df$qvalues <- NULL
	# df$rank <- NULL
	# df$leading_edge <- NULL
	# df$core_enrichment <- NULL
	col_order <- c("X", "p.adjust", "NES")
	df <- df[, col_order]
	df <- df %>% remove_rownames %>% column_to_rownames(var="X")
	colnames(df) <- c(vs_name, 'NES')
	df <- df[df[,vs_name] < 0.05, , drop = FALSE]
	
	# We plug in positive or negative log(padj) depending on the direction 
	# of the genes
	# padj 0.002 -> log10 -> -3
	# so for up we do -(-3) = 3
	# for down the opposite
	# if (df$NES >= 0) {
	# 	df[,vs_name] <- -log10(df[, vs_name]) 
	# } else {
	# 	df[,vs_name] <- log10(df[, vs_name]) 
	# }
	df[,vs_name] <- ifelse(df$NES >= 0, -log10(df[, vs_name]), log10(df[, vs_name]))
	df$NES <- NULL
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

go_rbinded <- rbind(go_merged_up, go_merged_down)

write.table(go_rbinded, file = go, col.names = TRUE, row.names = TRUE,
            sep = "\t", quote = FALSE)

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

merged <- rbind(go_rbinded, gsea_merged)
merged[is.na(merged)] <- 0

pheatmap(merged, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, 
		 filename = plot, width = 15, height = 8)

write.table(merged, file = go_gsea, col.names = TRUE, row.names = TRUE,
            sep = "\t", quote = FALSE)

# merged_subset <- merged[merged$LMO_LMH > 7 & merged$LMX_LMH > 7,]

save.image("paperino2.Rdata")

subs <- subset(merged, merged$`LMO-LMH` >= 20.0 | merged$`LMX-LMH` >= 20.0 | merged$`LMO-LMX` >= 20.0)
subs2 <- subset(merged, merged$`LMO-LMH` <= -20 | merged$`LMX-LMH` <= -20 | merged$`LMO-LMX` <=-20)

subset <- rbind(subs, subs2)

# merged_rownames <- merged
# merged_rownames <- tibble::rownames_to_column(merged_rownames, "GO")


pheatmap(subset, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE,
         fontsize_col = 10, width = 15, height = 8, cellwidth = 10, file = plot_20)

