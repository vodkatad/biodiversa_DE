library(tidyverse)

meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f$sample_id_new <- gsub('-','.', meda_f$sample_id_R, fixed=TRUE)
meda_f <- filter(meda_f, grepl("LMH", type))

vsd <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'SYMBOL'
vsd <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_new, "SYMBOL")]

write.table(vsd, file =  "LMH_gene_genealogyall.tsv", quote = FALSE, sep = "\t", col.names = TRUE)
