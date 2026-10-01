library(dplyr)

ncbi <- snakemake@input[["ncbi_f"]]
ivan <- snakemake@input[["gtf_f"]]
counts_f <- snakemake@input[["counts"]]
res <- snakemake@output[["file_counts"]]
gene_len_f <- snakemake@input[["gene_len"]]
res_gene <- snakemake@output[["file_gene_len"]]

#ncbi <- "/home/mferri/symbol_synonyms_official.tsv"
ncbi <- read.table(
  ncbi,
  sep = "\t",
  header = TRUE,
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

#ivan <- "/home/mferri/gene_symbols_v27.txt"
ivan <- read.table(
  ivan,
  col.names = "V1",
  stringsAsFactors = FALSE
)

## --- build synonym -> official symbol table ---
syn <- do.call(rbind, lapply(seq_len(nrow(ncbi)), function(i) {
  if (is.na(ncbi$Synonyms[i]) || ncbi$Synonyms[i] %in% c("", "-")) {
    return(NULL)
  }
  data.frame(
    gene = unlist(strsplit(ncbi$Synonyms[i], "\\|")),
    symbol_ncbi = ncbi$Symbol[i],
    stringsAsFactors = FALSE
  )
})) %>%
  filter(gene %in% ivan$V1) %>%
  distinct()

## flag ambiguous synonyms: stesso alias -> più di un official symbol diverso
## (teniamo traccia, ma ora scegliamo comunque il primo match)
syn_ambig <- syn %>%
  distinct(gene, symbol_ncbi) %>%
  group_by(gene) %>%
  filter(n() > 1) %>%
  ungroup()

ambiguous_synonyms <- unique(syn_ambig$gene)

## --- classify each symbol from our GTF ---
ivan_ncbi <- data.frame(
  gene = ivan$V1,
  tipo = "non trovato",
  symbol_ncbi = NA_character_,
  stringsAsFactors = FALSE
)

if (any(duplicated(ivan_ncbi$gene))) {
  warning(sum(duplicated(ivan_ncbi$gene)), " gene symbol duplicati nel file GTF di input")
}

idx_official <- match(ivan_ncbi$gene, ncbi$Symbol)
ivan_ncbi$tipo[!is.na(idx_official)] <- "ufficiale"
ivan_ncbi$symbol_ncbi[!is.na(idx_official)] <-
  ncbi$Symbol[idx_official[!is.na(idx_official)]]

## per i sinonimi, prendi il primo match (comportamento di default di match())
idx_synonym <- match(ivan_ncbi$gene, syn$gene)
idx <- is.na(idx_official) & !is.na(idx_synonym)
ivan_ncbi$tipo[idx] <- "sinonimo"
ivan_ncbi$symbol_ncbi[idx] <- syn$symbol_ncbi[idx_synonym[idx]]

## marca (solo per reportistica) quali sinonimi erano ambigui,
## ma NON tocchiamo symbol_ncbi: resta il primo match già assegnato
idx_ambig <- ivan_ncbi$tipo == "sinonimo" & ivan_ncbi$gene %in% ambiguous_synonyms
ivan_ncbi$tipo[idx_ambig] <- "sinonimo_ambiguo_primo_match"

table(ivan_ncbi$tipo)

sinonimi <- ivan_ncbi %>% filter(tipo %in% c("sinonimo", "sinonimo_ambiguo_primo_match"))
sinonimi_ambigui <- ivan_ncbi %>% filter(tipo == "sinonimo_ambiguo_primo_match")

## per reportistica: tutti i candidati scartati per ciascun ambiguo
candidati_ambigui <- syn_ambig %>%
  filter(gene %in% sinonimi_ambigui$gene) %>%
  arrange(gene)

## --- load counts ---
#counts <- ""
counts <- read.table(counts_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts$Geneid <- gsub("H_", "", counts$Geneid)
rename_map <- setNames(sinonimi$symbol_ncbi, sinonimi$gene)

idx_rename <- counts$Geneid %in% names(rename_map)
cat("Geneid nel counts che saranno rinominati:", sum(idx_rename), "\n")

idx_rename_ambig <- counts$Geneid %in% sinonimi_ambigui$gene
cat("...di cui rinominati scegliendo il primo match tra più candidati ambigui:", sum(idx_rename_ambig), "\n")
if (sum(idx_rename_ambig) > 0) print(unique(counts$Geneid[idx_rename_ambig]))

## applica la rinomina
counts_fix <- counts
counts_fix$Geneid[idx_rename] <- rename_map[counts_fix$Geneid[idx_rename]]

## aggrega i duplicati creati dalla rinomina (sia i "molti ivan -> uno ncbi" fisiologici,
## sia le eventuali coincidenze accidentali dei primo-match ambigui), sommando i conteggi
counts_fix <- counts_fix %>%
  group_by(Geneid) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  as.data.frame()

cat("Righe counts originale:", nrow(counts), "\n")
cat("Righe counts dopo rinomina + aggregazione:", nrow(counts_fix), "\n")
cat("Differenza (= numero di fusioni avvenute):", nrow(counts) - nrow(counts_fix), "\n")

gz <- gzfile(res, "w")
write.table(counts_fix, file = gz, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
close(gz)
# fare la stessa cosa anche per il gene_len
#gene_len <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/gene_len"
gene_len <- read.table(gene_len_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_len$Geneid <- gsub("H_", "", gene_len$Geneid)

rename_map <- setNames(sinonimi$symbol_ncbi, sinonimi$gene)

idx_rename <- gene_len$Geneid %in% names(rename_map)
cat("Geneid nel counts che saranno rinominati:", sum(idx_rename), "\n")

idx_rename_ambig <- gene_len$Geneid %in% sinonimi_ambigui$gene
cat("...di cui rinominati scegliendo il primo match tra più candidati ambigui:", sum(idx_rename_ambig), "\n")
if (sum(idx_rename_ambig) > 0) print(unique(gene_len$Geneid[idx_rename_ambig]))

## applica la rinomina
gene_len_fix <- gene_len
gene_len_fix$Geneid[idx_rename] <- rename_map[gene_len_fix$Geneid[idx_rename]]

## aggrega i duplicati creati dalla rinomina (sia i "molti ivan -> uno ncbi" fisiologici,
## sia le eventuali coincidenze accidentali dei primo-match ambigui), sommando i conteggi
gene_len_fix <- gene_len_fix %>%
  group_by(Geneid) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  as.data.frame()

cat("Righe counts originale:", nrow(gene_len), "\n")
cat("Righe counts dopo rinomina + aggregazione:", nrow(gene_len_fix), "\n")
cat("Differenza (= numero di fusioni avvenute):", nrow(gene_len) - nrow(gene_len_fix), "\n")

write.table(gene_len_fix, file=res_gene, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
