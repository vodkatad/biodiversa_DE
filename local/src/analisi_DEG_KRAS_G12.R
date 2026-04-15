## analisi deg KRAS

## CRC1598

## MRTX + Cet vs NT

d <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC1598/M_C_vs_NT/CRC1598_treat_cutoff0.05-M_C.vs.NT.deseq2.tsv"
d <- read.table(d, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

up <- d %>% filter(padj < 0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- d %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < - 0.5849625)

## CRC0031 

## MRTX vs NT

mnt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/MRTX_vs_NT/CRC0031_treat_cutoff0.05-MRTX.vs.NT.deseq2.tsv"
mnt <- read.table(mnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

upmnt <- mnt %>% filter(padj < 0.05)
upmnt <- upmnt %>% filter(log2FoldChange > 0.5849625)
length(rownames(upmnt))

downmnt <- mnt %>% filter(padj < 0.05)
downmnt <- downmnt %>% filter(log2FoldChange < - 0.5849625)
length(rownames(downmnt))

## Cet + Tram vs NT

ctnt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/C_T_vs_NT/CRC0031_treat_cutoff0.05-C_T.vs.NT.deseq2.tsv"
ctnt <- read.table(ctnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ctnt$genes <- rownames(ctnt)

upctnt <- ctnt %>% filter(padj < 0.05)
upctnt <- upctnt %>% filter(log2FoldChange > 0.5849625)
length(rownames(upctnt))

downctnt <- ctnt %>% filter(padj < 0.05)
downctnt <- downctnt %>% filter(log2FoldChange < - 0.5849625)
length(rownames(downctnt))

ctntdeg <- ctnt %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5849625)
ctntdeg$genes <- rownames(ctntdeg)

## MRTX + Cet vs NT

mcnt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_vs_NT/CRC0031_treat_cutoff0.05-M_C.vs.NT.deseq2.tsv"
mcnt <- read.table(mcnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mcnt$genes <- rownames(mcnt)

upmcnt <- mcnt %>% filter(padj < 0.05)
upmcnt <- upmcnt %>% filter(log2FoldChange > 0.5849625)
length(rownames(upmcnt))

downmcnt <- mcnt %>% filter(padj < 0.05)
downmcnt <- downmcnt %>% filter(log2FoldChange < - 0.5849625)
length(rownames(downmcnt))

mcntdeg <- mcnt %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5849625)
mcntdeg$genes <- rownames(mcntdeg)

## MRTX + Cet + Tram vs NT

mctnt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_T_vs_NT/CRC0031_treat_cutoff0.05-M_C_T.vs.NT.deseq2.tsv"
mctnt <- read.table(mctnt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mctnt$genes <- rownames(mctnt)

upmctnt <- mctnt %>% filter(padj < 0.05)
upmctnt <- upmctnt %>% filter(log2FoldChange > 0.5849625)
length(rownames(upmctnt))

downmctnt <- mctnt %>% filter(padj < 0.05)
downmctnt <- downmctnt %>% filter(log2FoldChange < - 0.5849625)
length(rownames(downmctnt))

mctntdeg <- mctnt %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5849625)
mctntdeg$genes <- rownames(mctntdeg)

## MRTX + Cet + Alb vs NT

mcant <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/results/CRC0031/M_C_B_vs_NT/CRC0031_treat_cutoff0.05-M_C_B.vs.NT.deseq2.tsv"
mcant <- read.table(mcant, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

upmcant <- mcant %>% filter(padj < 0.05)
upmcant <- upmcant %>% filter(log2FoldChange > 0.5849625)
length(rownames(upmcant))

downmcant <- mcant %>% filter(padj < 0.05)
downmcant <- downmcant %>% filter(log2FoldChange < - 0.5849625)
length(rownames(downmcant))

mcantdeg <- mcant %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5849625)

## VENN

## valutazione a monte con intersezioni. 
## mcant, mctnt e mcnt hanno stranamente tutti gli stessi geni analizzati: 20545.
## faccio delle intersezioni per vedere quanti di questi 20545 sono condivisi.

intmcnt_mctnt <- intersect(rownames(mcnt), rownames(mctnt))
length(intmcnt_mctnt)

intmcnt_mcant <- intersect(rownames(mcnt), rownames(mcant))
length(intmcnt_mcant)

## Sono gli stessi identici geni mossi - il trattamento farmacologico con doppia o tripla interessa gli stessi 20545 geni
## è la stessa cosa pratically?

## intersezione con l'altra doppia

intctnt_mctnt <- intersect(rownames(ctnt), rownames(mctnt))
length(intctnt_mctnt)

intctnt_mcant <- intersect(rownames(ctnt), rownames(mcant))
length(intctnt_mcant)

doppia <- intersect(rownames(ctnt), rownames(mcnt))
length(doppia)

library(ggVennDiagram)

## UP

venn_list <- list(
  MRTX_Cet_up = rownames(upmcnt),
  MRTX_Cet_Tram_up = rownames(upmctnt),
  MRTX_Cet_Alb_up = rownames(upmcant)
)

ggVennDiagram(venn_list)

venn_list_up_cet_cetalb <- list(
  MRTX_Cet_up = rownames(upmcnt),
  MRTX_Cet_Alb_up = rownames(upmcant)
)

p <- ggVennDiagram(venn_list_up_cet_cetalb)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

venn_list_up_cet_cettram <- list(
  MRTX_Cet_up = rownames(upmcnt),
  MRTX_Cet_Tram_up = rownames(upmctnt)
)

p <- ggVennDiagram(venn_list_up_cet_cettram)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

## DOWN

venn_list_down <- list(
  MRTX_Cet_down = rownames(downmcnt),
  MRTX_Cet_Tram_down = rownames(downmctnt),
  MRTX_Cet_Alb_down = rownames(downmcant)
)

ggVennDiagram(venn_list_down)

venn_list_down_cet_cetalb <- list(
  MRTX_Cet_down = rownames(downmcnt),
  MRTX_Cet_Alb_down = rownames(downmcant)
)

p <- ggVennDiagram(venn_list_down_cet_cetalb)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

venn_list_down_cet_cettram <- list(
  MRTX_Cet_down = rownames(downmcnt),
  MRTX_Cet_Tram_down = rownames(downmctnt)
)

p <- ggVennDiagram(venn_list_down_cet_cettram)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

## ALL

venn_alb <- list(
  MRTX_Cet = rownames(mcntdeg),
  MRTX_Cet_Alb = rownames(mcantdeg)
)

p <- ggVennDiagram(venn_alb)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

venn_tram <- list(
  MRTX_Cet = rownames(mcntdeg),
  MRTX_Cet_Tram = rownames(mctntdeg)
)

p <- ggVennDiagram(venn_tram)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))

## DOPPIE

venn_doppie <- list(
  Cet_Tram = rownames(ctntdeg),
  MRTX_Cet = rownames(mcntdeg)
)

p <- ggVennDiagram(venn_doppie)
p + coord_cartesian(clip = "off") +
  theme(plot.margin = margin(20, 60, 20, 60))


## check NA mono

geni_na <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/KRAS_G12/CRC0031_geni_NA_mono_MRTX.tsv"
geni_na <- read.table(geni_na, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geni_na <- rownames(geni_na)

mcntna <- mcnt %>% filter(rownames(mcnt) %in% geni_na)
mctntna <- mctnt %>% filter(rownames(mctnt) %in% geni_na)
mcantna <- mcant %>% filter(rownames(mcant) %in% geni_na)

## scatters 

m_cet_m_cet_tram <- merge(mcntdeg, mctntdeg, by="genes")
names(m_cet_m_cet_tram)[names(m_cet_m_cet_tram)=="log2FoldChange.x"] <- "log2FoldChange_MRTX+CET.vs.NT" 
names(m_cet_m_cet_tram)[names(m_cet_m_cet_tram)=="log2FoldChange.y"] <- "log2FoldChange_MRTX+CET+TRAM.vs.NT" 
ggplot(m_cet_m_cet_tram, aes(`log2FoldChange_MRTX+CET.vs.NT`, `log2FoldChange_MRTX+CET+TRAM.vs.NT`))+geom_point()+geom_smooth(method=lm)

cet_tram_m_cet_tram <- merge(ctntdeg, mctntdeg, by="genes")
names(cet_tram_m_cet_tram)[names(cet_tram_m_cet_tram)=="log2FoldChange.x"] <- "log2FoldChange_CET+TRAM.vs.NT"
names(cet_tram_m_cet_tram)[names(cet_tram_m_cet_tram)=="log2FoldChange.y"] <- "log2FoldChange_MRTX+CET+TRAM.vs.NT"
ggplot(cet_tram_m_cet_tram, aes(`log2FoldChange_CET+TRAM.vs.NT`, `log2FoldChange_MRTX+CET+TRAM.vs.NT`))+geom_point()+geom_smooth(method = lm)

m_cet_m_cet_tram_inversi <- m_cet_m_cet_tram
m_cet_m_cet_tram_inversi$segno <- m_cet_m_cet_tram_inversi$`log2FoldChange_MRTX+CET.vs.NT`*m_cet_m_cet_tram_inversi$`log2FoldChange_MRTX+CET+TRAM.vs.NT`
m_cet_m_cet_tram_inversi <- m_cet_m_cet_tram_inversi %>% filter(segno < -0.5)

m_cet_m_cet_tram_all<- merge(mcnt, mctnt, by = "genes", all = TRUE)
names(m_cet_m_cet_tram_all)[names(m_cet_m_cet_tram_all)=="log2FoldChange.x"] <- "log2FoldChange_MRTX+CET.vs.NT" 
names(m_cet_m_cet_tram_all)[names(m_cet_m_cet_tram_all)=="log2FoldChange.y"] <- "log2FoldChange_MRTX+CET+TRAM.vs.NT" 

names(m_cet_m_cet_tram_all)[names(m_cet_m_cet_tram_all)=="padj.x"] <- "padj_MRTX+CET.vs.NT" 
names(m_cet_m_cet_tram_all)[names(m_cet_m_cet_tram_all)=="padj.y"] <- "padj_MRTX+CET+TRAM.vs.NT" 

m_cet_m_cet_tram_all$DEG_MRTX_CET <- ifelse(m_cet_m_cet_tram_all$`padj_MRTX+CET.vs.NT` < 0.05 & abs(m_cet_m_cet_tram_all$`log2FoldChange_MRTX+CET.vs.NT`) > 0.5849625, "DEG_MRTX_CET", "NO_DEG_MRT_CET")
m_cet_m_cet_tram_all$DEG_MRTX_CET_TRAM <- ifelse(m_cet_m_cet_tram_all$`padj_MRTX+CET+TRAM.vs.NT` < 0.05 & abs(m_cet_m_cet_tram_all$`log2FoldChange_MRTX+CET+TRAM.vs.NT`) > 0.5849625, "DEG_MRTX_CET_TRAM", "NO_DEG_MRT_CET_TRAM")


m_cet_m_cet_tram_all <- m_cet_m_cet_tram_all %>%
  mutate(DEG_combined = case_when(
    DEG_MRTX_CET == "DEG_MRTX_CET" &
      DEG_MRTX_CET_TRAM == "DEG_MRTX_CET_TRAM" ~ "both",
    
    DEG_MRTX_CET == "DEG_MRTX_CET" &
      DEG_MRTX_CET_TRAM != "DEG_MRTX_CET_TRAM" ~ "MRTX_CET_only",
    
    DEG_MRTX_CET != "DEG_MRTX_CET" &
      DEG_MRTX_CET_TRAM == "DEG_MRTX_CET_TRAM" ~ "MRTX_CET_TRAM_only",
    
    TRUE ~ "none"
  ))

ggplot(m_cet_m_cet_tram_all,
       aes(`log2FoldChange_MRTX+CET.vs.NT`,
           `log2FoldChange_MRTX+CET+TRAM.vs.NT`,
           color = DEG_combined)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    none = "grey70",
    MRTX_CET_only = "#1f77b4",
    MRTX_CET_TRAM_only = "#d62728",
    both = "#2ca02c"
  )) +
  theme_minimal()


cet_tram_m_cet_tram_all<- merge(ctnt, mctnt, by = "genes", all = TRUE)
names(cet_tram_m_cet_tram_all)[names(cet_tram_m_cet_tram_all)=="log2FoldChange.x"] <- "log2FoldChange_CET+TRAM.vs.NT" 
names(cet_tram_m_cet_tram_all)[names(cet_tram_m_cet_tram_all)=="log2FoldChange.y"] <- "log2FoldChange_MRTX+CET+TRAM.vs.NT" 

names(cet_tram_m_cet_tram_all)[names(cet_tram_m_cet_tram_all)=="padj.x"] <- "padj_CET+TRAM.vs.NT" 
names(cet_tram_m_cet_tram_all)[names(cet_tram_m_cet_tram_all)=="padj.y"] <- "padj_MRTX+CET+TRAM.vs.NT" 

cet_tram_m_cet_tram_all$DEG_CET_TRAM <- ifelse(cet_tram_m_cet_tram_all$`padj_CET+TRAM.vs.NT` < 0.05 & abs(cet_tram_m_cet_tram_all$`log2FoldChange_CET+TRAM.vs.NT`) > 0.5849625, "DEG_CET_TRAM", "NO_DEG_CET_TRAM")
cet_tram_m_cet_tram_all$DEG_MRTX_CET_TRAM <- ifelse(cet_tram_m_cet_tram_all$`padj_MRTX+CET+TRAM.vs.NT` < 0.05 & abs(cet_tram_m_cet_tram_all$`log2FoldChange_MRTX+CET+TRAM.vs.NT`) > 0.5849625, "DEG_MRTX_CET_TRAM", "NO_DEG_MRT_CET_TRAM")

cet_tram_m_cet_tram_all <- cet_tram_m_cet_tram_all %>%
  mutate(DEG_combined = case_when(
    DEG_CET_TRAM == "DEG_CET_TRAM" &
      DEG_MRTX_CET_TRAM == "DEG_MRTX_CET_TRAM" ~ "both",
    
    DEG_CET_TRAM == "DEG_CET_TRAM" &
      DEG_MRTX_CET_TRAM != "DEG_MRTX_CET_TRAM" ~ "DEG_CET_TRAM_only",
    
    DEG_CET_TRAM != "DEG_CET_TRAM" &
      DEG_MRTX_CET_TRAM == "DEG_MRTX_CET_TRAM" ~ "MRTX_CET_TRAM_only",
    
    TRUE ~ "none"
  ))

ggplot(cet_tram_m_cet_tram_all,
       aes(`log2FoldChange_CET+TRAM.vs.NT`,
           `log2FoldChange_MRTX+CET+TRAM.vs.NT`,
           color = DEG_combined)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    none = "grey70",
    DEG_CET_TRAM_only = "#1f77b4",
    MRTX_CET_TRAM_only = "#d62728",
    both = "#2ca02c"
  )) +
  theme_minimal()
