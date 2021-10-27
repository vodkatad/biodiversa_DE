### heatmap shutte

library(pheatmap)
library(tidyverse)

# we have 6 GO and 3 GSEA
list_enrich_results <- snakemake@input
plot <- snakemake@output[["hm"]]

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

pheatmap(merged, cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 4.7,
         fontsize_col = 6, filename = plot, width = 15, height = 8)

# write.table merged

# merged_subset <- merged[merged$LMO_LMH > 7 & merged$LMX_LMH > 7,]




# go_Lmo_lmx <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE.tsv"
# gsea_lmo_lmx <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GSEA_results_H_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE.tsv"

# go_ox <- read.table(go_Lmo_lmx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# gs_ox <- read.table(gsea_lmo_lmx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# go_lmo_lmh <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMH.tsv"
# gsea_lmo_lmh <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GSEA_results_H_type_cutoff0.05-LMO_BASALE.vs.LMH.tsv"

# go_oh <- read.table(go_lmo_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# gs_oh <- read.table(gsea_lmo_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# go_lmx_lmh <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GO_results_type_cutoff0.05-LMX_BASALE.vs.LMH.tsv"
# gsea_lmx_lmh <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/GSEA_results_H_type_cutoff0.05-LMX_BASALE.vs.LMH.tsv"

# go_xh <- read.table(go_lmx_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# gs_xh <- read.table(gsea_lmx_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# # GO

# go_ox$X <- paste(go_ox$ID, go_ox$Description)
# go_ox$ID <- NULL
# go_ox$Description <- NULL
# go_ox$GeneRatio <- NULL
# go_ox$BgRatio <- NULL
# go_ox$pvalue <- NULL
# go_ox$qvalue <- NULL
# go_ox$geneID <- NULL
# go_ox$Count <- NULL
# col_order <- c("X", "p.adjust")
# go_ox <- go_ox[, col_order]
# go_ox <- go_ox %>% remove_rownames %>% column_to_rownames(var="X")
# colnames(go_ox) <- c("p.adjust_LMO_LMX")

# go_oh$X <- paste(go_oh$ID, go_oh$Description)
# go_oh$ID <- NULL
# go_oh$Description <- NULL
# go_oh$GeneRatio <- NULL
# go_oh$BgRatio <- NULL
# go_oh$pvalue <- NULL
# go_oh$qvalue <- NULL
# go_oh$geneID <- NULL
# go_oh$Count <- NULL
# col_order <- c("X", "p.adjust")
# go_oh <- go_oh[, col_order]
# go_oh <- go_oh %>% remove_rownames %>% column_to_rownames(var="X")
# colnames(go_oh) <- c("p.adjust_LMO_LMH")

# go_xh$X <- paste(go_xh$ID, go_xh$Description)
# go_xh$ID <- NULL
# go_xh$Description <- NULL
# go_xh$GeneRatio <- NULL
# go_xh$BgRatio <- NULL
# go_xh$pvalue <- NULL
# go_xh$qvalue <- NULL
# go_xh$geneID <- NULL
# go_xh$Count <- NULL
# col_order <- c("X", "p.adjust")
# go_xh <- go_xh[, col_order]
# go_xh <- go_xh %>% remove_rownames %>% column_to_rownames(var="X")
# colnames(go_xh) <- c("p.adjust_LMX_LMH")

# #GSEA
# gs_ox$ID <- NULL
# gs_ox$Description <- NULL
# gs_ox$setSize <- NULL
# gs_ox$enrichmentScore <- NULL
# gs_ox$NES <- NULL
# gs_ox$pvalue <- NULL
# gs_ox$qvalues <- NULL
# gs_ox$rank <- NULL
# gs_ox$leading_edge <- NULL
# gs_ox$core_enrichment <- NULL
# colnames(gs_ox) <- c("p.adjust_LMO_LMX")

# gs_oh$ID <- NULL
# gs_oh$Description <- NULL
# gs_oh$setSize <- NULL
# gs_oh$enrichmentScore <- NULL
# gs_oh$NES <- NULL
# gs_oh$pvalue <- NULL
# gs_oh$qvalues <- NULL
# gs_oh$rank <- NULL
# gs_oh$leading_edge <- NULL
# gs_oh$core_enrichment <- NULL
# colnames(gs_oh) <- c("p.adjust_LMO_LMH")

# gs_xh$ID <- NULL
# gs_xh$Description <- NULL
# gs_xh$setSize <- NULL
# gs_xh$enrichmentScore <- NULL
# gs_xh$NES <- NULL
# gs_xh$pvalue <- NULL
# gs_xh$qvalues <- NULL
# gs_xh$rank <- NULL
# gs_xh$leading_edge <- NULL
# gs_xh$core_enrichment <- NULL
# colnames(gs_xh) <- c("p.adjust_LMX_LMH")

# go_gs_ox <- rbind(go_ox, gs_ox)
# go_gs_oh <- rbind(go_oh, gs_oh)
# go_gs_xh <- rbind(go_xh, gs_xh)

# go_gs_ox <- go_gs_ox[go_gs_ox$p.adjust_LMO_LMX < 0.05, , drop = FALSE]
# go_gs_ox$p.adjust_LMO_LMX <- -log10(go_gs_ox$p.adjust_LMO_LMX) 
# go_gs_oh <- go_gs_oh[go_gs_oh$p.adjust_LMO_LMH < 0.05, , drop = FALSE]
# go_gs_oh$p.adjust_LMO_LMH <- -log10(go_gs_oh$p.adjust_LMO_LMH)
# go_gs_xh <- go_gs_xh[go_gs_xh$p.adjust_LMX_LMH < 0.05, , drop = FALSE]
# go_gs_xh$p.adjust_LMX_LMH <- -log10(go_gs_xh$p.adjust_LMX_LMH)

# merged_ox_oh <- merge(go_gs_ox, go_gs_oh, by = "row.names", all = TRUE)
# merged_ox_oh <- merged_ox_oh %>% remove_rownames %>% column_to_rownames(var="Row.names")

# merged <- merge(merged_ox_oh, go_gs_xh, by = "row.names", all = TRUE)
# merged <- merged %>% remove_rownames %>% column_to_rownames(var="Row.names")

# merged[is.na(merged)] <- 0

# pheatmap(merged, cluster_rows = TRUE, cluster_cols = FALSE)
