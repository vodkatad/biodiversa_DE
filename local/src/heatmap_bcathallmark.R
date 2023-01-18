sample <- "/scratch/trcanmed/DE_RNASeq/dataset/deg_tcf7l2_mut/samples_data"
sample <- read.table(sample, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi <- rownames(sample)

genes <- read.xlsx("/scratch/trcanmed/DE_RNASeq/dataset/deg_tcf7l2_mut/hallmark_wnt_bcat.xlsx")

vsd <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz", quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'

genisign <- genes$Member  

vsd_subset <- vsd %>% filter(symbol %in% genisign)
rownames(vsd_subset) <- vsd_subset$symbol
vsd_subset$symbol <- NULL
vsd_subset <- as.data.frame(t(vsd_subset))
vsd_subset$sample <- rownames(vsd_subset)
vsd_subset <- vsd_subset %>% filter(sample %in% casi)
vsd_subset$sample <- NULL

pheatmap(vsd_subset)
