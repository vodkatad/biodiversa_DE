## analysis signatures wt e mut

wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_vsd.tsv.gz"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt$genes <- rownames(wt)

i <- "/home/mferri/GO0060337_type_I_interferon-mediated_signaling_pathway.tsv"
i <- read.table(i, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
i <- unique(i$SYMBOL)

#GOBP_TYPE_II_INTERFERON_PRODUCTION from GSEA
ii2 <- c("MIR708","KLRC4-KLRK1","EBI3","CD96","CEBPG","CD226","LILRB1",
         "ARID5A","LILRB4","RIPK3","BTN3A2","BTN3A1","CD160","PGLYRP2","PGLYRP3","SLAMF6",
         "CCR7","CR1","SIRPA","IL23R","ZFPM1","DDIT3","NLRP6","F2RL1","KLRK1","SCRIB","IL27",
         "ABL1","ZNF683","PTPN22","GAS6","GATA3","IL36RN","IFNL1","LGALS9B","HMSD","PYCARD",
         "CD274","HLA-A","HLA-DPA1","HLA-DPB1","HLA-DRB1","HMGB1","HRAS","HSPD1","IRF8","IRGM",
         "APP","IL1B","IL1R1","IL2","IL10","IL12A","IL12B","IL12RB1","IL12RB2","IL18","INHA",
         "INHBA","ISL1","JAK2","LGALS9","LTA","MIR24-1","FOXP3","TLR7","TLR8","PDE4B","PDE4D",
         "UFC1","IL23A","CYRIB","CD244","IL20RB","TLR9","SASH3","UFSP2","AXL","PRNP","CRTAM",
         "HMHB1","IL21","RARA","TRIM27","BCL3","XCL1","VSIR","CLEC7A","SLAMF1","LGALS9C",
         "SLC11A1","C1QBP","TLR3","TLR4","TNF","TNFSF4","CCR2","TXK","TYK2","SCGB1A1","WNT5A",
         "ZP3","LAPTM5","FZD5","UBA5","ZC3H12A","PDCD1LG2","CD276","SLC7A5","HAVCR2","RIPK2","FADD",
         "IL18R1","PGLYRP1","IL33","CD2","CD3E","IL1RL1","CD14","IL27RA","CD47","ISG15","NR1H4")

found_in_go <- c("RIPK3", "VSIR","F2RL1","LGALS9","PDE4D","ISG15","IL23R","IL33","CD244","LGALS9C","IL36RN","IL18R1","IL1B","NR1H4","LGALS9B","PGLYRP3","SIRPA","TNF"
                 ,"HLA-DRB1","SLC11A1","NLRP6","TLR3","INHBA","IFNL1","IL1RL1","EBI3")

check <- intersect(ii2, found_in_go)
intersect_real <- intersect(i, ii2)

wt_i <- wt %>% filter(genes %in% i)
wt_i$genes <- NULL
wt_i$mean <- rowMeans(wt_i)

deg <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg_i <- deg %>% filter(genes %in% i)

go_wt_up <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
go_wt_up <- read.table(go_wt_up, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_wt_up <- go_wt_up %>% filter(p.adjust < 0.05)

go_mut_up <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv"
go_mut_up <- read.table(go_mut_up, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_mut_up <- go_mut_up %>% filter(p.adjust < 0.05)

common <- intersect(go_wt_up$ID, go_mut_up$ID)

for (i in rownames(go_wt_up)) {
  if (i %in% common) {
    go_wt_up[i,"common"] <- "yes"
  } else {
    go_wt_up[i,"common"] <- "no"
  }
}

for (i in rownames(go_mut_up)) {
  if (i %in% common) {
    go_mut_up[i, "common"] <- "yes"
  } else {
    go_mut_up[i, "common"] <- "no"
  }
}

go_wt_up <- go_wt_up[order(go_wt_up$p.adjust),]
go_mut_up <- go_mut_up[order(go_mut_up$p.adjust),]

write.xlsx(go_wt_up, file = "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_up_N2.vs.NE_significative_confronto_MUT.xlsx")
write.xlsx(go_mut_up, file = "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_up_N2.vs.NE_significative_confronto_WT.xlsx")

gas <- "/home/mferri/GO0140896_cGAS:STING_signaling_pathway.tsv"
gas <- read.table(gas, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gas <- unique(gas$SYMBOL)

deg_magno <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
deg_magno <- read.table(deg_magno, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg_magno$genes <- rownames(deg_magno)
deg_gas <- deg_magno %>% filter(genes %in% gas)
deg_gas <- deg_gas %>% filter(padj < 0.05)
deg_gas <- deg_gas %>% filter(abs(log2FoldChange) > 0.5849625)

ad <- "/home/mferri/GO_0007155_cell_adhesion.tsv"
ad <- read.table(ad, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ad <- unique(ad$SYMBOL)

apical <- "/home/mferri/GO_0045177_apical_part_of_cell.tsv"
apical <- read.table(apical, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
apical <- unique(apical$SYMBOL)

intersect(ad, apical)


## CHECK GENES ALL WNT

## check all WNT signatures in GO

wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt <- wt %>% filter(p.adjust < 0.05)
wt <- wt[grepl("Wnt", wt$Description), ]
geni_wt <- as.list(wt$geneID)
geni_wt <- unlist(geni_wt)
geni_wt <- strsplit(geni_wt, split = "/")
geni_wt <- unique(unlist(geni_wt))

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mut <- mut %>% filter(p.adjust < 0.05)
mut <- mut[grepl("Wnt", mut$Description), ]
geni_mut <- as.list(mut$geneID)
geni_mut <- unlist(geni_mut)
geni_mut <- strsplit(geni_mut, split = "/")
geni_mut <- unique(unlist(geni_mut))

all <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv"
all <- read.table(all, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
all <- all %>% filter(p.adjust < 0.05)
all <- all[grepl("Wnt", all$Description), ]
geni_all <- as.list(all$geneID)
geni_all <- unlist(geni_all)
geni_all <- strsplit(geni_all, split = "/")
geni_all <- unique(unlist(geni_all))

setdiff(geni_wt, geni_mut)
# [1] "PTK7"     "LGR5"     "LGR6"     "TMEM131L" "DISC1"    "ADGRA2"   "WNT10B"   "APOE"     "NID1"    
# [10] "FZD10"    "FZD2"     "DLX5"     "DACT3"    "GPRC5B"   "LRP4"     "TPBGL"    "DKKL1"    "TMEM88B" 
# [19] "ARHGEF19"
setdiff(geni_mut, geni_wt)
# [1] "FZD9"  "RSPO4" "BARX1"
setdiff(geni_all, geni_wt)
#[1] "BARX1" -> gli all sono guidati dai wt
setdiff(geni_wt, geni_all)
# [1] "TMEM131L" "DISC1"    "WNT10B"   "FZD2"     "GPRC5B"   "LRP4"     "TPBGL"    "DKKL1"    "TMEM88B" 
# [10] "ARHGEF19"
setdiff(geni_all, geni_mut)
# [1] "PTK7"   "LGR5"   "LGR6"   "APOE"   "NID1"   "ADGRA2" "DLX5"   "DACT3"  "FZD10" 
setdiff(geni_mut, geni_all)
# [1] "FZD9"  "RSPO4"

geni <- unique(c(geni_all, geni_mut, geni_wt))

lfc_wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
lfc_wt <- read.table(lfc_wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lfc_wt$genes <- rownames(lfc_wt)
names(lfc_wt)[names(lfc_wt)=="log2FoldChange"] <- "log2FoldChange_wt"
lfc_wt <- lfc_wt %>% filter(genes %in% geni)
lfc_wt <- lfc_wt[,c("genes", "log2FoldChange_wt")]

lfc_mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
lfc_mut <- read.table(lfc_mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lfc_mut$genes <- rownames(lfc_mut)
names(lfc_mut)[names(lfc_mut)=="log2FoldChange"] <- "log2FoldChange_mut"
lfc_mut <- lfc_mut %>% filter(genes %in% geni)
lfc_mut <- lfc_mut[,c("genes", "log2FoldChange_mut")]

merged <- merge(lfc_wt, lfc_mut, by="genes")
rownames(merged) <- merged$genes
for (i in rownames(merged)) {
  gene_values <- data.frame(
    Condizione = c("WT", "MUT"),
    log2FoldChange = c(merged[i, "log2FoldChange_wt"], merged[i, "log2FoldChange_mut"])
  )
  p <- ggplot(gene_values, aes(x = Condizione, y = log2FoldChange)) +
    geom_boxplot() + geom_jitter(height = 0)
  
  print(p)
}

## CHECK citochine globale
cit <- "/home/mferri/GO_0019221_cytokine-mediated_signaling_pathway.tsv"
cit <- read.table(cit, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cit <- unique(cit$SYMBOL)

deg <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/geno_cutoff0.05-N2.vs.NE.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg <- deg %>% filter(genes %in% cit)
deg <- deg %>% filter(padj < 0.05)
deg <- deg %>% filter(abs(log2FoldChange) > 0.5849625)

write.table(deg[,"genes"], file="/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/citochine", quote = FALSE,
            sep = "\t", col.names = FALSE, row.names = FALSE)

## CHECK GENI SIGNATURE PER BOXPLOT SSGSEA
jakstat <-  
