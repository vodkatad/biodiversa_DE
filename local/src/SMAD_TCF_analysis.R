## smad tcf

## mutazioni
df <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/smads/smads_specific-mutations_list_GOI_PDX_alltiers.tsv"
df <- read.table(df, quote = "",sep = "\t", header = TRUE, stringsAsFactors = FALSE)

uni <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/smads/smads_specific-mutations_list_PDX_alltiers.tsv"
uni <- read.table(uni, quote = "",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
uni <- uni %>% filter(!mutated_gene %in% c("SMAD3", "SMAD4", "TCF7L2"))
uni$type <- "WT"
uni <- uni[,c("smodel", "type")]
uni <- uni[!duplicated(uni$smodel),]

length(intersect(df$smodel, uni$smodel))
uni <- uni %>% filter(!smodel %in% df$smodel)

smad_models <- unique(df$smodel[df$mutated_gene %in% c("SMAD3", "SMAD4")])
tcf_models  <- unique(df$smodel[df$mutated_gene == "TCF7L2"])

mut_status <- df %>%
  select(smodel, mutated_gene) %>%
  distinct() %>%
  mutate(mutated = 1) %>%
  pivot_wider(names_from = mutated_gene, values_from = mutated, values_fill = 0)

mut_status <- mut_status %>%
  mutate(SMAD_altered = ifelse((`SMAD3` == 1 | `SMAD4` == 1), 1, 0),
         SMAD_group = ifelse(SMAD_altered == 1, "SMAD_mut", "SMAD_wt"))

uni$type <- NULL
uni$SMAD4 <- 0
uni$TCF7L2 <- 0
uni$SMAD3 <- 0 
uni$SMAD_altered <- 0
uni$SMAD_group <- "SMAD_wt"

mut_status <- rbind(mut_status, uni)

table(mut_status$SMAD_group, mut_status$TCF7L2)

mut_status %>%
  group_by(SMAD_group) %>%
  summarise(freq_TCF7L2 = mean(TCF7L2),
            n_models = n())
tbl <- table(mut_status$SMAD_group, mut_status$TCF7L2)
#16% smad mutati nei wt di tcf
#28% smad mutati nei mut di tcf
fisher.test(tbl)

order_dep_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx" 
order_dep <- read_xlsx(order_dep_f) 
order_dep <- order_dep %>% filter(!case == "CRC0152")

merged <- order_dep %>%
  rename(smodel = case) %>%            
  merge(mut_status, by = "smodel")   

merged %>%
  group_by(SMAD_group) %>%
  summarise(
    mean_dep = mean(co_comp_N2, na.rm = TRUE),
    median_dep = median(co_comp_N2, na.rm = TRUE),
    n = n()
  )

ggplot(merged, aes(x = SMAD_group, y = co_comp_N2, fill = SMAD_group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2)

wilcox.test(co_comp_N2 ~ SMAD_group, data = merged)

## delezioni

xeno <- "/mnt/cold1/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_thresholded.by_genes.txt"
xeno <- read.table(xeno, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## dubbio su quei casi che hanno smad3 1 e smad4 -1 come considerarli es CRC0081
classify_cn <- function(x) {
  if (x <= -1 ) return("loss")
  if (x >= 1)  return("gain")
  return("neutral")
}

smad_genes <- xeno %>%
  filter(Gene.Symbol %in% c("SMAD4","SMAD3")) %>%
  select(-Locus.ID, -Cytoband)

smad_cn <- smad_genes[,-1]

smad_class <- smad_cn
for(i in 1:nrow(smad_cn)){
  for(j in 1:ncol(smad_cn)){
    smad_class[i,j] <- classify_cn(smad_cn[i,j])
  }
}

smad_status <- character(ncol(smad_class))
for(j in 1:ncol(smad_class)){
  col_vals <- smad_class[,j]
  if ("loss" %in% col_vals) {
    smad_status[j] <- "loss"
  } else if ("gain" %in% col_vals) {
    smad_status[j] <- "gain"
  } else {
    smad_status[j] <- "neutral"
  }
}

tcf_gene <- xeno %>%
  filter(Gene.Symbol == "TCF7L2") %>%
  select(-Locus.ID, -Cytoband)

tcf_cn <- tcf_gene[,-1]

tcf_status <- character(ncol(tcf_cn))
for(j in 1:ncol(tcf_cn)){
  tcf_status[j] <- classify_cn(tcf_cn[1,j])
}

samples_xeno <- colnames(xeno)[-(1:3)]
sample_short <- str_extract(samples_xeno, "CRC\\d+")

cn_table <- data.frame(
  case = sample_short,
  SMAD_status = smad_status,
  TCF7L2_status = tcf_status
)

merged_cn <- merge(cn_table, order_dep, by="case")

merged_cn$SMAD_group <- ifelse(merged_cn$SMAD_status == "loss", 
                               "SMAD_loss", "SMAD_neutral")
merged_cn$TCF7L2_group <- ifelse(merged_cn$SMAD_status == "loss", 
                               "SMAD_loss", "SMAD_neutral")
wilcox.test(co_comp_N2 ~ SMAD_group, data = merged_cn)

ggplot(merged_cn, aes(x = SMAD_group, y = co_comp_N2)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.7)) +
  geom_jitter(width = 0.1, size = 2) 

## espressione
## fare lo stesso plot con la divisione wt e mut per tcf

expr <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/fpkm.tsv.gz"
expr <- read.table(expr, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
expr$genes <- rownames(expr)
expr <- expr %>% filter(genes %in% c("SMAD3", "SMAD4"))
expr$genes <- NULL
col_groups <- sub("_R[0-9]+$", "", colnames(expr))

groups <- split(seq_along(colnames(expr)), col_groups)

expr_mean <- sapply(groups, function(idx) {
  if(length(idx) == 1) {
    expr[, idx]
  } else {
    apply(expr[, idx, drop=FALSE], 1, mean)
  }
})

expr_mean <- as.data.frame(expr_mean)

expr_NE <- expr_mean[, grepl("_NE$", colnames(expr_mean))]
colnames(expr_NE) <- gsub("_NE", "", colnames(expr_NE))
expr_NE <- as.data.frame(t(expr_NE))
expr_NE$SMAD_mean <- rowMeans(expr_NE[, c("SMAD3", "SMAD4")], na.rm = TRUE)
expr_NE$case <- rownames(expr_NE)

merged_expr <- merge(expr_NE, order_dep, by="case")

cor.test(merged_expr$SMAD_mean, merged_expr$co_comp_N2)

ggplot(merged_expr, aes(x = SMAD_mean, y = co_comp_N2)) +
  geom_point() +
  geom_smooth(method="lm")

mut <- c("CRC0148", "CRC1331", "CRC0399", "CRC1278", "CRC0277", 
         "CRC0327", "CRC1729", "CRC0152", "CRC0196", "CRC0059", 
         "CRC0065", "CRC0464", "CRC0316", "CRC1239") 

rownames(merged_expr) <- merged$smodel

for (i in rownames(merged_expr)) {
  if (i %in% mut) {
    merged_expr[i,"type"] <- "MUT"
  } else {
    merged_expr[i, "type"] <- "WT"
  }
}

ggplot(merged_expr, aes(x = SMAD_mean, y = co_comp_N2)) +
  geom_point(aes(color=type)) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = case))

wt <- merged_expr[merged_expr$type == "WT",]
cor.test(wt$SMAD_mean, wt$co_comp_N2)

mut <- merged_expr[merged_expr$type == "MUT",]
cor.test(mut$SMAD_mean, mut$co_comp_N2)

## aggiungere alla tabella delle mutazioni anche chi ha una delezione
merged_all <- merge(mut_status, cn_table[, c("case", "SMAD_status")],
                by.x = "smodel", by.y = "case", all.x = TRUE)
merged_all <- merged_all[!duplicated(merged_all$smodel),]
rownames(merged_all) <- merged_all$smodel

merged_all$SMAD_status_all <- ifelse(
  merged_all$SMAD_group == "SMAD_mut" | merged_all$SMAD_status %in% c("loss","gain"),
  "MUT",
  "WT"
)

tbl <- table(merged_all$SMAD_status_all, merged_all$TCF7L2)
tbl
fisher.test(tbl)
