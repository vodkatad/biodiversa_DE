## arricchimenti mutazionali 

d <- "/home/mferri/Domain_TCF.csv"
d <- read.csv(d)
d <- d[,c("Description","Start","End")]
domain <- d$Description
d <- as.data.frame(d)
togliere <- c("HMG-box_TCF7-like", 
              "Catenin binding domain superfamily", "hmgende2",
              "HMG_box", "HMG_BOX_2",
              "HMG-box", "TRANSCRIPTION FACTOR 7 FAMILY MEMBER")
d <- d %>% dplyr::filter(!Description %in% togliere)
d$Description <- as.character(d$Description)
d$Description <- c("C-Clamp", "CTNNB1_binding", "HMG-box")

m <- read.table("/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/TCF7L2_SOX9/manual/mutinfo_added_back.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(m) <- c("genomic", "sample_id", "source", "annotation")

# dopo osservazioni su 
#chr10:113165541:TCCCCTGTTTCTAGGAGAAA:T	CRC1182LMX0A01001TUMD08000V2	WES_PriMets	splice_acceptor_variant&coding_sequence_variant&intron_variant:ENST00000627217.3:c.1392-13_1397del	NA	splice_acceptor_variant&coding_sequence_variant&intron_variant	TRUE
# può diventare p.r464spl
idx <- which(m$sample_id == "CRC1182LMX0A01001TUMD08000V2")
m$annotation[idx] <- "splice_acceptor_variant&coding_sequence_variant&intron_variant:ENST00000627217.3:p.R464spl"
m <- m %>%
  mutate(
    aa_pos = ifelse(
      str_detect(annotation, "p\\.[A-Z][0-9]+"),
      as.numeric(str_extract(annotation, "(?<=p\\.[A-Z])([0-9]+)")),
      NA
    ),
    mut_type = str_extract(annotation, "^[^:]+")
  )

m$domain <- "Inter-domain"

for(i in 1:nrow(d)) {
  sel <- !is.na(m$aa_pos) & m$aa_pos >= d$Start[i] & m$aa_pos <= d$End[i]
  m$domain[sel] <- d$Description[i]
}

m <- m[,c("sample_id", "domain")]
m$model <- substr(m$sample_id, 1, 7)
m <- m[!duplicated(m$model),]

uni <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/TCF7L2_SOX9/universe_models.tsv"
uni <- read.table(uni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
uni <- uni[!duplicated(uni$smodel),, drop=FALSE]
bb <- "/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/badboys"
bb <- read.table(bb)
uni <- uni %>% filter(!smodel %in% bb$V1)
wt <- uni
wt$type <- "WT"
wt <- wt %>% filter(!smodel %in% m$model)

muttcf <- m
names(muttcf)[names(muttcf)=="model"] <- "smodel"
muttcf$type <- "MUT"
muttcfdomain <- muttcf
muttcf <- muttcf[,c("smodel", "type")]
muttcf <- muttcf[!duplicated(muttcf$smodel),]
tcf <- rbind(muttcf, wt)
muttcfdomain <- muttcfdomain[,c("smodel", "domain")]
names(muttcfdomain)[names(muttcfdomain)=="domain"] <- "type"
muttcfdomain <- rbind(muttcfdomain, wt)
wt_models <- unique(wt$smodel)
mut_models <- unique(m$model)
all_models <- unique(c(wt_models, mut_models))

n_wt <- length(wt_models)
n_mut_total <- length(mut_models)

domains <- unique(m$domain)

results_list <- list()

for (d in domains) {
  # mutati in questo dominio
  mut_in_d <- unique(m$model[m$domain == d])
  n_in_d <- length(mut_in_d)
  # mutati in altri domini
  n_mut_not_in_d <- n_mut_total - n_in_d
  tbl <- matrix(c(n_in_d, n_mut_not_in_d,
                  0,      n_wt),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(
                  Status = c("MUT", "WT"),
                  Domain = c("In_domain", "Not_in_domain")
                ))
  ft <- fisher.test(tbl)
  results_list[[d]] <- data.frame(
    domain = d,
    n_mut_in_domain = n_in_d,
    n_mut_not_in_domain = n_mut_not_in_d,
    n_wt = n_wt,
    p.value = ft$p.value,
    odds_ratio = if(!is.null(ft$estimate)) ft$estimate else NA
  )
}

results_df <- do.call(rbind, results_list)

## arricchimenti mutazioni varie
top10 <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/top10_3/top10_3/top10_3_specific-mutations_list_GOI_PDX_wtiers.tsv"
top10 <- read.table(top10, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bb <- "/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/badboys"
bb <- read.table(bb)
top10 <- top10 %>% filter(!smodel %in% bb$V1)

df_list <- split(top10, top10$mutated_gene)

df_list <- lapply(df_list, function(x) {
  out <- unique(x["smodel"])
  out$type <- "MUT"
  out
})

uni <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/mutcheck/ba_round/TCF7L2_SOX9/universe_models.tsv"
uni <- read.table(uni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
uni <- uni[!duplicated(uni$smodel),, drop=FALSE]
uni <- uni %>% filter(!smodel %in% bb$V1)

df_complete <- lapply(names(df_list), function(g) {
  mut_df <- df_list[[g]]
  wt_df <- uni[!(uni$smodel %in% mut_df$smodel), , drop = FALSE]
  wt_df$type <- "WT"
  out <- rbind(mut_df, wt_df)
  out$gene <- g   # aggiungi colonna con il gene di riferimento
  rownames(out) <- NULL
  out
})

names(df_complete) <- names(df_list)

fisher_results <- list()

for (g in names(df_complete)) {
  gene_df <- df_complete[[g]]
  mut_gene <- gene_df$smodel[gene_df$type == "MUT"]
  wt_gene  <- gene_df$smodel[gene_df$type == "WT"]
  mut_tcf <- tcf$smodel[tcf$type == "MUT"]
  wt_tcf  <- tcf$smodel[tcf$type == "WT"]
  a <- length(intersect(wt_tcf, wt_gene))   # WT_TCF & WT_gene
  b <- length(intersect(wt_tcf, mut_gene))  # WT_TCF & MUT_gene
  c <- length(intersect(mut_tcf, wt_gene))  # MUT_TCF & WT_gene
  d <- length(intersect(mut_tcf, mut_gene)) # MUT_TCF & MUT_gene
  tbl <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(
                  TCF = c("WT_TCF", "MUT_TCF"),
                  Gene = c("WT_gene", "MUT_gene")
                ))
  ft <- fisher.test(tbl)
  fisher_results[[g]] <- data.frame(
    gene = g,
    a = a, b = b, c = c, d = d,
    p.value = ft$p.value,
    odds_ratio = if (!is.null(ft$estimate)) ft$estimate else NA
  )
}

fisher_df <- do.call(rbind, fisher_results)
rownames(fisher_df) <- NULL

genivstcf <- fisher_df

## fisher tra dominio e altre mutazioni

domains <- setdiff(unique(muttcfdomain$type), "WT")
results <- list()

for (d in domains) {
  mut_domain <- muttcfdomain$smodel[muttcfdomain$type == d]
  wt_domain  <- muttcfdomain$smodel[muttcfdomain$type == "WT"]
  for (g in names(df_complete)) {
    gene_df <- df_complete[[g]]
    mut_gene <- gene_df$smodel[gene_df$type == "MUT"]
    wt_gene  <- gene_df$smodel[gene_df$type == "WT"]
    a <- length(intersect(wt_domain, wt_gene))   # WT_domain & WT_gene
    b <- length(intersect(wt_domain, mut_gene))  # WT_domain & MUT_gene
    c <- length(intersect(mut_domain, wt_gene))  # MUT_domain & WT_gene
    dval <- length(intersect(mut_domain, mut_gene)) # MUT_domain & MUT_gene
    tbl <- matrix(c(a, b, c, dval), nrow = 2,
                  dimnames = list(
                    Domain = c("WT_domain", "MUT_domain"),
                    Gene   = c("WT_gene", "MUT_gene")
                  ))
    ft <- fisher.test(tbl)
    results[[paste(g, d, sep = "_")]] <- data.frame(
      gene = g,
      domain = d,
      a = a, b = b, c = c, d = dval,
      p.value = ft$p.value,
      odds_ratio = if (!is.null(ft$estimate)) ft$estimate else NA
    )
  }
}

fisher_results <- do.call(rbind, results)
rownames(fisher_results) <- NULL

genivsdomini <- fisher_results

## risposta ba su arid1a

a <- top10
a <- top10 %>% filter(mutated_gene == "ARID1A")
a <- a[!duplicated(a$smodel),]
a$tipo <- ifelse(grepl("missense", a$mut_type), "missenso", "LoF")
a <- a[,c("smodel","tipo")]
merged <- merge(a, tcf, by="smodel", all.x=TRUE, all.y=TRUE)
merged$tipo[is.na(merged$tipo)] <- "WT"
colnames(merged) <- c("smodel", "ARID1A", "TCF7L2")
merged$ARID1A <- gsub("LoF", "MUT", merged$ARID1A)
merged$ARID1A <- gsub("missenso", "MUT", merged$ARID1A)
rownames(merged) <- merged$smodel
for (i in rownames(merged)) {
  if (merged[i, "ARID1A"] == "MUT" & merged[i, "TCF7L2"]=="MUT") {
    merged[i,"mutation"] <- "both"
  } else if (merged[i, "ARID1A"] == "MUT" & merged[i, "TCF7L2"]=="WT") {
    merged[i,"mutation"] <- "ARID1A"
  } else if (merged[i, "ARID1A"] == "WT" & merged[i, "TCF7L2"]=="MUT") {
    merged[i,"mutation"] <- "TCF7L2"
  } else {
    merged[i,"mutation"] <- "WT"
  }
}

write.table(merged, file="ARID1A_TCF_mut.tsv", quote = FALSE, sep = "\t",col.names = TRUE, row.names = FALSE)

arid1a <- "/home/mferri/excel_geni_tcf_bcat/ARID1A_TCF_mut.tsv"
a <- read.table(arid1a, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tab <- dplyr::count(a, ARID1A, TCF7L2)
tab <- tab %>%
  mutate(
    class = case_when(
      ARID1A == "WT"  & TCF7L2 == "WT"  ~ "WT / WT",
      ARID1A == "MUT" & TCF7L2 == "WT"  ~ "ARID1A only",
      ARID1A == "WT"  & TCF7L2 == "MUT" ~ "TCF7L2 only",
      ARID1A == "MUT" & TCF7L2 == "MUT" ~ "Double mutant"
    )
  )
fill_cols <- c(
  "WT / WT"        = "#d9d9d9",  # grigio chiaro
  "ARID1A only"    = "#4c72b0",  # blu
  "TCF7L2 only"    = "#dd8452",  # arancio
  "Double mutant"  = "#55a868"   # verde
)
tab <- tab %>%
  mutate(
    ARID1A = factor(ARID1A, levels = c("WT", "MUT")),
    TCF7L2 = factor(TCF7L2, levels = c("WT", "MUT"))
  )
ggplot(tab, aes(x = ARID1A, y = TCF7L2, fill = class)) +
  geom_tile(color = "black", linewidth = 1) +
  geom_text(aes(label = n), size = 6, fontface = "bold") +
  scale_fill_manual(values = fill_cols) +
  scale_y_discrete(limits = rev(levels(tab$TCF7L2))) +
  labs(
    x = "ARID1A",
    y = "TCF7L2"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  )
