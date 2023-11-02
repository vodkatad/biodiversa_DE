## pre oncoprint e waterfall top 10 biobanca

library(tidyverse)

lmx_f <- snakemake@input[["merged"]]
res <- snakemake@output[["result"]]
vaf_tsv <- snakemake@output[["vaf_result"]]
prot_tsv <- snakemake@output[["prot_result"]]

#lmx_f <- "/scratch/trcanmed/snakegatk/dataset/biobanca_targeted_pdx/mutect/merged_longformat_wtiers.tsv"
lmx <- read.table(lmx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx <- lmx %>% filter(!af == 0.000)

mut_x <- table(lmx$gene)
mut_x <- sort(mut_x, decreasing = TRUE)
mut_x <- mut_x[1:10]
mut_x <- as.data.frame(mut_x)
mut_x <- mut_x$Var1
mut_x <- as.character(mut_x)
mut_x <- c(mut_x, "CTNNB1", "NRAS", "BRAF")

lmx <- read.table(lmx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx[c('type', 'def', "mutation", "esone", "mut_id")] <- str_split_fixed(lmx$cds, ':', 5)
lmx <- lmx %>% filter(gene %in% mut_x)

df <- lmx
df$mut_id <- NULL
df$cds <- NULL
df$type <- NULL
df$def <- NULL
df$mutation <- NULL
df$esone <- NULL
df <- df %>% filter(!af==0)
# for (i in seq(rownames(df))) {
#   if (df[i,"af"]==0.0) {
#     df[i, "af_mat"] <- 0
#   } else {
#     df[i, "af_mat"] <- 1
#   }
# }
df$model <- substr(df$lgenealogy, 1, 7)
df$af <- 1
df$lgenealogy <- NULL
df <- df[!duplicated(df), ]
#df$af <- NULL
result <- as.data.frame(df %>% pivot_wider(names_from = model, values_from = af))
result[is.na(result)] <- 0
rownames(result) <- result$gene
result$gene <- NULL
## aggiungere CRC1472 al fondo con tutti zero
result$CRC1472 <- 0
result <- as.data.frame(t(result))

write.table(result, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

vaf <- lmx
vaf <- as.data.frame(vaf[,c(2,4,5)])
vaf$model <- substr(vaf$lgenealogy, 1, 7)
vaf$lgenealogy <- NULL
vaf <- vaf %>% filter(!af==0)
vaf_res <- as.data.frame(vaf %>% pivot_wider(names_from = model, values_from = af))
vaf_res <- as.data.frame(lapply(vaf_res, gsub, pattern='c', replacement=''))
vaf_res <- as.data.frame(lapply(vaf_res, gsub, pattern='[()]', replacement=''))
rownames(vaf_res) <- vaf_res$gene
vaf_res$gene <- NULL
vaf_res <- data.frame(lapply(vaf_res, as.character), stringsAsFactors=FALSE, row.names = rownames(vaf_res))
vaf_res[vaf_res == 'NULL'] <- 0

write.table(vaf_res, file=vaf_tsv, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

prot <- lmx
prot <- prot[,c(1,2,4,5)]
prot <- prot %>% filter(!af==0)
prot$af <- NULL
prot$model <- substr(prot$lgenealogy, 1, 7)
prot$lgenealogy <- NULL
prot_res <- as.data.frame(prot %>% pivot_wider(names_from = model, values_from = mut_id))
prot_res <- as.data.frame(lapply(prot_res, gsub, pattern='c', replacement=''))
prot_res <- as.data.frame(lapply(prot_res, gsub, pattern='[()]', replacement=''))
prot_res <- as.data.frame(lapply(prot_res, gsub, pattern='"', replacement=''))
prot_res <- as.data.frame(lapply(prot_res, gsub, pattern='"', replacement=''))
rownames(prot_res) <- prot_res$gene
prot_res$gene <- NULL
prot_res <- data.frame(lapply(prot_res, as.character), stringsAsFactors=FALSE, row.names = rownames(prot_res))
prot_res[prot_res == 'NULL'] <- "WT"

write.table(prot_res, file=prot_tsv, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
