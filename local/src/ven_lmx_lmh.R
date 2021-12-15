## venn diagram 4434 (stromal) vs up x-h
library(readxl)
library(ggvenn)

#rifare questo venn coi down

## intersect btw universe and stromal
g <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/type_cutoff0.05-LMX_BASALE.vs.LMH.deseq2.tsv"
genes_uni <- read.table(g, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes<- rownames(genes_uni)
rownames(genes_uni) <- NULL
genes_uni <- cbind(genes_uni,genes)
genes_uni <- genes_uni %>% mutate(genes = gsub("H_", "", genes))
col_order <- c("genes", "log2FoldChange")
genes_uni<- genes_uni[, col_order]


g <- genes_uni %>% remove_rownames %>% column_to_rownames(var="genes")

universe <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/type_cutoff0.05-LMX_BASALE.vs.LMH.gouniverse.tsv"
univ <- read.table(universe, quote = "", sep = "\t", header = FALSE)
u <- univ$V1
stromal <- read_excel("/scratch/trcanmed/DE_RNASeq/local/share/data/Stromal_genes.xlsx")
s <- stromal$Stromal_Genes


int <- as.data.frame(intersect(u, s))
write.table(int, file = "inteserct_stromal_universe.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
sgint <- int$`intersect(u, s)`

pca_100 <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/100_geni_pca.tsv"
pca <- read.table(pca_100, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

x <- list("Universo_Claudio" = sgint, "Pca100" = pca$V1)

ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)
setwd("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/")
ggsave("venn_diagram_xh.png")
genes_uni_lfc <- genes_uni[genes_uni$genes %in% int$`intersect(u, s)`,]



# meda <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
# meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# meda_f <- filter(meda_f, grepl("LMO_BASALE", type))
# 
# vsd <- "/scratch/trcanmed//DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
# vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
# genes <- rownames(vsd)
# rownames(vsd) <- NULL
# vsd <- cbind(genes,vsd)
# vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
# colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'
# vsd_subset <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_R, "symbol")]
# 
# vsd_subset <- vsd_subset %>% remove_rownames %>% column_to_rownames(var="symbol")
# 
# sds <- apply(vsd_subset, 1, sd)
# zerosd <- (names(sds[sds==0]))
# colnames(zerosd) <- "GENES"
# vsd_subset <- tibble::rownames_to_column(vsd_subset, "GENES")
# 
# intnopdo <- intersect(s, zerosd)

## venn diagram with up x-h

downxh <- "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/type_cutoff0.05-LMX_BASALE.vs.LMH.goinsplit_down.tsv"
dxh <- read.table(downxh, quote = "", sep = "\t", header = FALSE )
xh <- dxh$V1

genes_uxh <- genes_uni[genes_uni$genes %in% dxh$V1,]

setwd("/home/mferri/")
plot(density(genes_uni_lfc$log2FoldChange), main="Density plot genes universe")
ggsave("Density_plot_genes_universe.png")

# png("Density_plot_genes_universe.png", width=1280, height=720, units="px", type="cairo")
# plot(density(genes_uni_lfc$log2FoldChange))
# graphics.off()
# 
setwd("/home/mferri/")

png("Density_plot_genes_down_XH.png", width=1280, height=720, units="px", type="cairo")
plot(density(genes_uxh$log2FoldChange))
graphics.off()

x <- list("Universo_Claudio" = sgint, "Down_X_H" = xh)

ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)
setwd("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Cutoff0.05_LFC0.584/")
ggsave("venn_diagram_xh.png")

png("venn_diagram_xh.png", width=1280, height=720, units="px", type="cairo")
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 4)
graphics.off()


stromal_genes <- intersect(xh, sgint)
stromal_genes_f <- data.frame(matrix(nrow = length(stromal_genes), ncol = 1))
for (i in seq(1, length(stromal_genes))){
  stromal_genes_f[i,1] <- paste0("H_", stromal_genes[i])
}
names(stromal_genes_f)[names(stromal_genes_f) == "matrix.nrow...length.stromal_genes...ncol...1."] <- 'genes'
write.table(stromal_genes_f, file = "stromal_genes.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

phyper(84, 31788 - 2893, 2896, 1589, lower.tail = FALSE)
