## check expr fingerprinting

lmo <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
meda <- "/mnt/cold1/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"

lmo <- "/home/mferri/Ulisse_cold1/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
meda <- "/home/mferri/Ulisse_cold1/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"

meda <- read.table(meda, sep='\t', quote="", header=TRUE)
meda$sample_id <- gsub('-', '.', meda$sample_id, fixed = TRUE)
meda <- meda %>%
  filter(type %in% c("LMO_BASALE", "LMO_BASALE.1", "LMO_BASALE.2"))

lmo <- read.table(lmo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(lmo) <- str_remove(rownames(lmo), 'H_')
lmo <- as.data.frame(t(lmo))
lmo$genealogy <- rownames(lmo)
lmo <- lmo %>% filter(genealogy %in% meda$sample_id)
lmo$samples <- substr(rownames(lmo), 1, 7)

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/vsd.tsv.gz"
vsd <- "/home/mferri/Ulisse_cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/vsd.tsv.gz"
vsd <- "//home/mferri/Ulisse_cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/vsd.tsv.gz"
vsd <- read.table(vsd, sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)
#names(vsd)[names(vsd)=="CRC0542_NE_R2_"] <- "CRC0542_NE_R2"
colnames(vsd) <- gsub("sh_", "", colnames(vsd))
vsd_mean <- vsd %>%
  as_tibble() %>%
  mutate(gene = rownames(vsd)) %>%
  relocate(gene) %>%
  pivot_longer(
    -gene,
    names_to = "sample",
    values_to = "value"
  ) %>%
  mutate(sample_base = str_remove(sample, "_R[0-9]+$")) %>%
  group_by(gene, sample_base) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sample_base,
    values_from = mean_value
  )

vsd_mean <- as.data.frame(vsd_mean)
rownames(vsd_mean) <- vsd_mean$gene
vsd_mean$gene <- NULL

samples_ne <- substr(colnames(vsd_mean), 1, 7)
lmo <- lmo %>% filter(samples %in% samples_ne)
lmo$genealogy <- NULL

rn <- lmo$samples

replicate_index <- ave(seq_along(rn), rn, FUN = seq_along)
new_names <- paste0(rn, "_R", replicate_index)
rownames(lmo) <- new_names
lmo$samples <- NULL
lmo <- as.data.frame(t(lmo))

lmo_mean <- lmo %>%
  as_tibble() %>%
  mutate(gene = rownames(lmo)) %>%
  relocate(gene) %>%
  pivot_longer(
    -gene,
    names_to = "sample",
    values_to = "value"
  ) %>%
  mutate(sample_base = str_remove(sample, "_R[0-9]+$")) %>%
  group_by(gene, sample_base) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sample_base,
    values_from = mean_value
  )

lmo_mean <- as.data.frame(lmo_mean)
rownames(lmo_mean) <- lmo_mean$gene
lmo_mean$gene <- NULL

## tolgo CRC0456 perché non c'è in biobanca
#vsd_mean$CRC0456_NE <- NULL
colnames(lmo_mean) <- paste0(colnames(lmo_mean), "_BASALE")

bcat <- vsd_mean
bcat <- as.data.frame(t(bcat))
bcat$type <- substr(rownames(bcat), 9,10)
bcat <- bcat %>% filter(type %in% c("sc"))
bcat$type <- NULL
bcat <- as.data.frame(t(bcat))
vsd_mean <- bcat

vsd_mean$genes <- rownames(vsd_mean)
lmo_mean$genes <- rownames(lmo_mean)

d <- merge(vsd_mean, lmo_mean, by="genes")
rownames(d) <- d$genes
d$genes <- NULL

### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
means <- apply(d, 1, mean)
med <- median(means)

de <- d[means > med,]

sds <- apply(de, 1, sd)
### there are some very high sd?
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### 0.1086  0.3540  0.4702  0.5313  0.6411  3.7110 
# noisy_genes_thr <- quantile(sds, 0.9) # tenere o no? non � una differenza cos� enorme come prima...toglilo

# de <- de[sds < noisy_genes_thr,]
# sds <- sds[sds < noisy_genes_thr]

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

desd1 <- desd
names(desd1) <- substr(names(desd1), 1, 14)

lmo <- desd1[grepl('BASALE', names(desd1),, fixed=TRUE)]
bcat <- desd1[grepl('scr', names(desd1),, fixed=TRUE)]

colnames(lmo) <- substr(colnames(lmo),0,7)
colnames(bcat) <- substr(colnames(bcat),0,7)
clmo <- lmo[, names(bcat)]
cbcat <- bcat[, names(clmo)]
### check also for genes
all_genes <- intersect(rownames(lmo), rownames(bcat))
clmo <- clmo[all_genes,]
cbcat <- cbcat[all_genes,]


if (all(colnames(clmo)!=colnames(cbcat)) & all(rownames(clmo)!=rownames(cbcat))) {
  stop('Brutto llama!')
}

res <- cor(clmo, cbcat)
colnames(res) <- paste0(colnames(res), "_BASALE")
rownames(res) <- paste0(rownames(res), "_scr")

pheatmap(res, cluster_rows = FALSE, cluster_cols = FALSE)

#write.table(res, gzfile(opt$output), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

biob <- "/home/mferri/Ulisse_cold1/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMX-LMO_correlation_simo_buoni.tsv.gz"
biob <- read.table(biob, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
biob <- as.matrix(biob)

diag_bcat <- diag(res)
diag_bio <- diag(biob)
