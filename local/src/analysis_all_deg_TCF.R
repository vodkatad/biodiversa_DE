## Analisi all DEG

## N2.vs.NE general
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/table_cases.tsv",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg2 <- deg %>% filter(padj < 0.05)
deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
write.xlsx(deg2, file="geni_differenziali_N2.vs.NE.xlsx", rowNames=TRUE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv",
                 quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GO_up_N2.vs.NE_singificative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## MUT_N2.vs.NE
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg2 <- deg %>% filter(padj < 0.05)
deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
write.xlsx(deg2, file="geni_differenziali_MUT_N2.vs.NE.xlsx", rowNames=TRUE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_up_N2.vs.NE_significative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## WT_N2.vs.NE
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg2 <- deg %>% filter(padj < 0.05)
deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
write.xlsx(deg2, file="geni_differenziali_WT_N2.vs.NE.xlsx", rowNames=TRUE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_up_N2.vs.NE_significative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_results_geno_cutoff0.05-N2.vs.NE_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## basali MUT vs WT
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/mut_samples_data_NE",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/geno_mut_cutoff0.05-NE_MUT.vs.NE_WT.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/GO_results_geno_mut_cutoff0.05-NE_MUT.vs.NE_WT_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
#go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/GO_results_geno_mut_cutoff0.05-NE_MUT.vs.NE_WT_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## sh bcat b2 vs scr
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/table_cases.tsv",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))
s <- s %>% filter(!geno == "t2")

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg2 <- deg %>% filter(padj < 0.05)
deg2 <- deg2 %>% filter(abs(log2FoldChange) > 0.5849625)
write.xlsx(deg2, file="geni_differenziali_bcat.vs.scr.xlsx", rowNames=TRUE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/GO_results_geno_cutoff0.05-b2.vs.scr_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/GO_up_bcat.vs.scr_significative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/GO_results_geno_cutoff0.05-b2.vs.scr_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(pvalue < 0.05)
go_down <- go_down[order(go_down$p.adjust),]
write.xlsx(go_down,"/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/GO_down_bcat.vs.scr_significative.xlsx")
go_down <- head(go_down, 20)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## check overlap mut.vs.wt con globale N2.vs.NE
deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main/geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up_glo <- deg %>% filter(padj<0.05)
up_glo <- up_glo %>% filter(log2FoldChange > 0.5849625)
down_glo <- deg %>% filter(padj < 0.05)
down_glo <- down_glo %>% filter(log2FoldChange < -0.5849625)

geno <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/geno_mut_cutoff0.05-NE_MUT.vs.NE_WT.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up_geno <- geno %>% filter(padj < 0.05)
up_geno <- up_geno %>% filter(log2FoldChange >0.5849625)

down_geno <- geno %>% filter(padj < 0.05)
down_geno <- down_geno %>% filter(log2FoldChange < - 0.5849625)

library(VennDiagram)
library(grid)

genes_up_glo <- rownames(up_glo)
genes_up_geno <- rownames(up_geno)
genes_down_glo <- rownames(down_glo)
genes_down_geno <- rownames(down_geno)

# Venn 1: up_glo vs up_geno
venn_up <- venn.diagram(
  x = list(Global_N2.vs.NE = genes_up_glo, Basali_MUT.vs.WT = genes_up_geno),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_up)

# Venn 2: down_glo vs down_geno
venn_down <- venn.diagram(
  x = list(Global_N2.vs.NE = genes_down_glo, Basali_MUT.vs.WT = genes_down_geno),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_down)


## check MUT N2.vs.NE vs WT_N2.vs.NE
mut <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv", 
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up_mut <- mut %>% filter(padj < 0.05)
up_mut <- up_mut %>% filter(log2FoldChange >0.5849625)

down_mut <- mut %>% filter(padj < 0.05)
down_mut <- down_mut %>% filter(log2FoldChange < - 0.5849625)

wt <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                 quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up_wt <- wt %>% filter(padj < 0.05)
up_wt <- up_wt %>% filter(log2FoldChange >0.5849625)

down_wt <- wt %>% filter(padj < 0.05)
down_wt <- down_wt %>% filter(log2FoldChange < - 0.5849625)

genes_up_mut <- rownames(up_mut)
genes_up_wt <- rownames(up_wt)
genes_down_mut <- rownames(down_mut)
genes_down_wt <- rownames(down_wt)


# Venn 1: up_mut vs up_wt
venn_up <- venn.diagram(
  x = list(MUT_N2.vs.NE = genes_up_mut, WT_N2.vs.NE = genes_up_wt),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_up)

# Venn 2: down_mut vs down_wt
venn_down <- venn.diagram(
  x = list(MUT_N2.vs.NE = genes_down_mut, WT_N2.vs.NE = genes_down_wt),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_down)


### Confronto significativi BCAT.vs.SCR vs significativi TCF small N2.vs.NE
bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE)
bcat <- bcat %>% filter(padj < 0.05)
bcat <- bcat %>% filter(abs(log2FoldChange) > 0.5849625)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
bcat$genes <- rownames(bcat)

tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main_with_samples_bcat/geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
for (i in rownames(tcf_small)) {
  if (abs(tcf_small[i,"log2FoldChange_tcf_small_N2.vs.NE"]) > 0.5849625 & tcf_small[i,"padj"] < 0.05) {
    tcf_small[i,"differenziali_tcf"] <- "YES"
  } else {
    tcf_small[i,"differenziali_tcf"] <- "NO"
  }
}
tcf_small$genes <- rownames(tcf_small)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(x=log2FoldChange_bcat.vs.scr, y=log2FoldChange_tcf_small_N2.vs.NE))+geom_point()+xlim(-5,5)+ylim(-5,5)

merged$direction_bcat <- ifelse(merged$log2FoldChange_bcat.vs.scr > 0, "up", "down")
merged$direction_tcf <- ifelse(merged$differenziali_tcf == "YES", "sing", "no")

tab <- table(merged$direction_bcat, merged$direction_tcf)

fisher.test(tab)

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_tcf")]

## boxplot
bcat$direction <- ifelse(bcat$log2FoldChange_bcat.vs.scr > 0, "up", "down")
bcat$type <- "BCAT"

tcf_small$direction <- ifelse(tcf_small$log2FoldChange_tcf_small_N2.vs.NE > 0, "up", "down")
tcf_small$type <- "TCF"

bcat_df <- bcat[, c("genes", "log2FoldChange_bcat.vs.scr", "direction", "type")]
colnames(bcat_df)[2] <- "log2FoldChange"

tcf_df <- tcf_small[, c("genes", "log2FoldChange_tcf_small_N2.vs.NE", "direction", "type")]
colnames(tcf_df)[2] <- "log2FoldChange"

combined <- rbind(bcat_df, tcf_df)

combined$group <- paste0(combined$type, "_", combined$direction)

combined$group <- factor(combined$group, levels = c("BCAT_up",  "TCF_up","BCAT_down", "TCF_down"))

ggplot(combined, aes(x = group, y = log2FoldChange, fill = direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("up" = "tomato", "down" = "steelblue")) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "log2 Fold Change", fill = "Direction") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


write.xlsx(merged, file="globale_significativi_bcat_differenziali_tcf.xlsx")

# differenziali TCF 
cor_yes <- cor.test(
  merged$log2FoldChange_bcat.vs.scr[merged$differenziali_tcf == "YES"],
  merged$log2FoldChange_tcf_small_N2.vs.NE[merged$differenziali_tcf == "YES"])

#  NON differenziali 
cor_no <- cor.test(
  merged$log2FoldChange_bcat.vs.scr[merged$differenziali_tcf == "NO"],
  merged$log2FoldChange_tcf_small_N2.vs.NE[merged$differenziali_tcf == "NO"])

merged_genes <- merged
rownames(merged_genes) <- merged_genes$genes
merged_genes <- merged_genes %>% filter(differenziali_tcf == "YES")
for (i in rownames(merged_genes)) {
  if (merged_genes[i, "log2FoldChange_bcat.vs.scr"] > 0 & merged_genes[i, "log2FoldChange_tcf_small_N2.vs.NE"]<0) {
    merged_genes[i, "keep"] <- "yes"
  } else if (merged_genes[i, "log2FoldChange_bcat.vs.scr"] < 0 & merged_genes[i, "log2FoldChange_tcf_small_N2.vs.NE"]>0) {
    merged_genes[i, "keep"] <- "yes"
  } else {
    merged_genes[i, "keep"] <- "no"
  }
}
merged_genes <- merged_genes %>% filter(keep == "yes")
merged_genes$keep <- NULL
merged_genes <- merged_genes[,c("log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE")]

merged_genes <- merged_genes %>% rownames_to_column("Genes")
colnames(merged_genes) <- c("Genes", "BCAT", "TCF")

df_long <- merged_genes %>%
  pivot_longer(cols = c("BCAT","TCF"), names_to = "Condition", values_to = "logFC")

ggplot(df_long, aes(x = Condition, y = logFC, color = Condition, group = Genes)) +
  geom_point(size = 3) +
  geom_line(aes(group = Genes), color = "gray60") +
  facet_wrap(~Genes, scales = "free_y") +
  theme_bw() +
  labs(y = "log2FC")

# vsd_bcat <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/vsd.tsv.gz"
# vsd_bcat <- read.table(vsd_bcat, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# vsd_bcat$genes <- rownames(vsd_bcat)
# vsd_bcat <- vsd_bcat %>% filter(genes %in% rownames(merged_genes))
# vsd_bcat$genes <- NULL
# vsd_bcat <- as.data.frame(t(vsd_bcat))
# vsd_bcat$t2 <- substr(rownames(vsd_bcat), 12, 13)
# vsd_bcat <- vsd_bcat %>% filter(!t2 == "t2")
# vsd_bcat$t2 <- NULL
# rownames(vsd_bcat) <- gsub("sh_", "", rownames(vsd_bcat))
# 
# vsd_tcf <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main_with_samples_bcat/vsd.tsv.gz"
# vsd_tcf <- read.table(vsd_tcf, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# vsd_tcf$genes <- rownames(vsd_tcf)
# vsd_tcf <- vsd_tcf %>% filter(genes %in% rownames(merged_genes))
# vsd_tcf$genes <- NULL
# vsd_tcf <- as.data.frame(t(vsd_tcf))



### Confronto significativi TCF small N2.vs.NE vs significativi BCAT.vs.SCR
tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_main_with_samples_bcat/geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tcf_small <- tcf_small %>% filter(padj < 0.05)
tcf_small <- tcf_small %>% filter(abs(log2FoldChange) > 0.5849625)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
tcf_small$genes <- rownames(tcf_small)

bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
for (i in rownames(bcat)) {
  if (abs(bcat[i,"log2FoldChange_bcat.vs.scr"]) > 0.5849625 & bcat[i,"padj"] < 0.05) {
    bcat[i,"differenziali_bcat"] <- "YES"
  } else {
    bcat[i,"differenziali_bcat"] <- "NO"
  }
}
bcat$genes <- rownames(bcat)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(y=log2FoldChange_bcat.vs.scr, x=log2FoldChange_tcf_small_N2.vs.NE))+geom_point()+xlim(-7,7)+ylim(-7,7)

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_bcat")]

write.xlsx(merged, file="globale_significativi_tcf_differenziali_bcat.xlsx")

# differenziali TCF 
cor_yes_tcf_bcat <- cor.test(
  merged$log2FoldChange_bcat.vs.scr[merged$differenziali_bcat == "YES"],
  merged$log2FoldChange_tcf_small_N2.vs.NE[merged$differenziali_bcat == "YES"])

#  NON differenziali 
cor_no_tcf_bcat <- cor.test(
  merged$log2FoldChange_bcat.vs.scr[merged$differenziali_bcat == "NO"],
  merged$log2FoldChange_tcf_small_N2.vs.NE[merged$differenziali_bcat == "NO"])


tcf_only <- merged %>% filter(differenziali_bcat == "NO")
tcf_only_up <- tcf_only %>% filter(log2FoldChange_tcf_small_N2.vs.NE > 0.5849625)
tcf_only_up$log2FoldChange_bcat.vs.scr <- NULL

write.xlsx(tcf_only_up, file="privati_tcf_up.xlsx")
tcf_only_down <- tcf_only %>% filter(log2FoldChange_tcf_small_N2.vs.NE < - 0.5849625)
tcf_only_down$log2FoldChange_bcat.vs.scr <- NULL
write.xlsx(tcf_only_down, file="privati_tcf_down.xlsx")

merged_genes <- merged
rownames(merged_genes) <- merged_genes$genes
merged_genes <- merged_genes %>% filter(differenziali_bcat == "YES")
for (i in rownames(merged_genes)) {
  if (merged_genes[i, "log2FoldChange_bcat.vs.scr"] > 3 & merged_genes[i, "log2FoldChange_tcf_small_N2.vs.NE"]>2) {
    merged_genes[i, "keep"] <- "yes"
  } else {
    merged_genes[i, "keep"] <- "no"
  }
}
merged_genes <- merged_genes %>% filter(keep == "yes")
merged_genes$keep <- NULL
merged_genes <- merged_genes[,c("log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE")]

## Confronto WT e MUT separati per bcat.vs.scr
wt_bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_bcat.vs.scr/WT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wt_bcat$genes <- rownames(wt_bcat)

up_wt <- wt_bcat %>% filter(padj < 0.05)
up_wt <- up_wt %>% filter(log2FoldChange >0.5849625)

down_wt <- wt_bcat %>% filter(padj < 0.05)
down_wt <- down_wt %>% filter(log2FoldChange < - 0.5849625)

mut_bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_bcat.vs.scr/MUT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                       quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mut_bcat$genes <- rownames(mut_bcat)
up_mut <- mut_bcat %>% filter(padj < 0.05)
up_mut <- up_mut %>% filter(log2FoldChange >0.5849625)

down_mut <- mut_bcat %>% filter(padj < 0.05)
down_mut <- down_mut %>% filter(log2FoldChange < - 0.5849625)

genes_up_mut <- rownames(up_mut)
genes_up_wt <- rownames(up_wt)
genes_down_mut <- rownames(down_mut)
genes_down_wt <- rownames(down_wt)


# Venn 1: up_mut vs up_wt
venn_up <- venn.diagram(
  x = list(MUT_BCAT.vs.SCR = genes_up_mut, WT_BCAT.vs.SCR = genes_up_wt),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_up)

# Venn 2: down_mut vs down_wt
venn_down <- venn.diagram(
  x = list(MUT_BCAT.vs.SCR = genes_down_mut, WT_BCAT.vs.SCR = genes_down_wt),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn_down)

## expr tcf7l2 basali

deg <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/geno_mut_cutoff0.05-NE_MUT.vs.NE_WT.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)
deg <- deg %>% filter(genes == "TCF7L2")

samples <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/mut_samples_data_NE"
samples <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_MUT.vs.WT/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$genes <- rownames(vsd)
vsd <- vsd %>% filter(genes == "TCF7L2")
vsd$genes <- NULL
vsd <- as.data.frame(t(vsd))
vsd$id <- rownames(vsd)

vsd <- merge(vsd, samples, by="id")
wt <- vsd %>% filter(geno_mut == "NE_WT")
mut <- vsd %>% filter(geno_mut == "NE_MUT")

wt_tcf <- mean(wt$TCF7L2)
mut_tcf <- mean(mut$TCF7L2)

## basali subset in biobanca
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/GO_results_geno_mut_cutoff0.05-MUT.vs.WT_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
#write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_up_N2.vs.NE_significative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF/GO_results_geno_mut_cutoff0.05-MUT.vs.WT_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))

## TCGA mut vs wt
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCGA/samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$Sample.ID))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCGA/geno_cutoff0.05-MUT.vs.WT.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

go_up <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCGA/GO_results_geno_cutoff0.05-MUT.vs.WT_up.tsv",
                    quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_up <- go_up %>% filter(p.adjust < 0.05)
go_up <- go_up[order(go_up$p.adjust),]
#write.xlsx(go_up, "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE/WT_GO_up_N2.vs.NE_significative.xlsx")
go_up <- head(go_up, 20)

df <- go_up
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


go_down <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCGA/GO_results_geno_cutoff0.05-MUT.vs.WT_down.tsv",
                      quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_down <- go_down %>% filter(p.adjust < 0.05)
go_down <- head(go_down, 20)

df <- go_down
x <- "GeneRatio"
df$num <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 1))
df$def <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[", 2))
df$x <- df$num/df$def
df$GeneRatio <- df$x

idx <- order(df[["GeneRatio"]], decreasing = TRUE)
df$Description <- factor(df$Description,
                         levels=rev(unique(df$Description[idx])))

ggplot(df, aes_string(x="GeneRatio", y="Description", size="Count", color="pvalue")) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = "pvalue",
                         guide=guide_colorbar(reverse=TRUE)) +
  ylab(NULL) +
  scale_size(range=c(3, 8))


## TCF but biobanca LMX
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF_LMX/samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF_LMX/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
up <- deg %>% filter(padj<0.05)
up <- up %>% filter(log2FoldChange > 0.5849625)

down <- deg %>% filter(padj < 0.05)
down <- down %>% filter(log2FoldChange < -0.5849625)

## extended biobanca
s <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_MUT.vs.WT_biobanca/samples_data",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
length(unique(s$model))

deg <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/basali_biobanca_like_TCF_LMX/geno_mut_cutoff0.05-MUT.vs.WT.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


### Confronto significativi BCAT.vs.SCR vs significativi TCF small N2.vs.NE WT
bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_bcat.vs.scr/WT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE)
bcat <- bcat %>% filter(padj < 0.05)
bcat <- bcat %>% filter(abs(log2FoldChange) > 0.5849625)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
bcat$genes <- rownames(bcat)

tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE_samples_bcat/WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
for (i in rownames(tcf_small)) {
  if (abs(tcf_small[i,"log2FoldChange_tcf_small_N2.vs.NE"]) > 0.5849625 & tcf_small[i,"padj"] < 0.05) {
    tcf_small[i,"differenziali_tcf"] <- "YES"
  } else {
    tcf_small[i,"differenziali_tcf"] <- "NO"
  }
}
tcf_small$genes <- rownames(tcf_small)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(x=log2FoldChange_bcat.vs.scr, y=log2FoldChange_tcf_small_N2.vs.NE, color=differenziali_tcf))+geom_point()+geom_smooth(method = lm)+ggtitle("WT")

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_tcf")]

write.xlsx(merged, file="WT_significativi_bcat_differenziali_tcf.xlsx")

### Confronto significativi TCF small N2.vs.NE vs significativi BCAT.vs.SCR WT
tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_N2.vs.NE_samples_bcat//WT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tcf_small <- tcf_small %>% filter(padj < 0.05)
tcf_small <- tcf_small %>% filter(abs(log2FoldChange) > 0.5849625)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
tcf_small$genes <- rownames(tcf_small)

bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_bcat.vs.scr/WT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
for (i in rownames(bcat)) {
  if (abs(bcat[i,"log2FoldChange_bcat.vs.scr"]) > 0.5849625 & bcat[i,"padj"] < 0.05) {
    bcat[i,"differenziali_bcat"] <- "YES"
  } else {
    bcat[i,"differenziali_bcat"] <- "NO"
  }
}
bcat$genes <- rownames(bcat)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(y=log2FoldChange_bcat.vs.scr, x=log2FoldChange_tcf_small_N2.vs.NE, color=differenziali_bcat))+geom_point()+geom_smooth(method = lm)+ggtitle("WT")

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_bcat")]

tcf_only_wt <- merged %>% filter(differenziali_bcat == "NO")
tcf_only_up_wt <- tcf_only_wt %>% filter(log2FoldChange_tcf_small_N2.vs.NE > 0.5849625)
tcf_only_up$log2FoldChange_bcat.vs.scr <- NULL

write.xlsx(tcf_only_up_wt, file="privati_tcf_wt_up.xlsx")
tcf_only_down_wt <- tcf_only_wt %>% filter(log2FoldChange_tcf_small_N2.vs.NE < - 0.5849625)
tcf_only_down_wt$log2FoldChange_bcat.vs.scr <- NULL
write.xlsx(tcf_only_down_wt, file="privati_tcf_wt_down.xlsx")


write.xlsx(merged, file="WT_significativi_tcf_differenziali_bcat.xlsx")

### Confronto significativi BCAT.vs.SCR vs significativi TCF small N2.vs.NE MUT
bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_bcat.vs.scr/MUT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE)
bcat <- bcat %>% filter(padj < 0.05)
bcat <- bcat %>% filter(abs(log2FoldChange) > 0.5849625)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
bcat$genes <- rownames(bcat)

tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE_samples_bcat/MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
for (i in rownames(tcf_small)) {
  if (abs(tcf_small[i,"log2FoldChange_tcf_small_N2.vs.NE"]) > 0.5849625 & tcf_small[i,"padj"] < 0.05) {
    tcf_small[i,"differenziali_tcf"] <- "YES"
  } else {
    tcf_small[i,"differenziali_tcf"] <- "NO"
  }
}
tcf_small$genes <- rownames(tcf_small)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(x=log2FoldChange_bcat.vs.scr, y=log2FoldChange_tcf_small_N2.vs.NE, color=differenziali_tcf))+geom_point()+geom_smooth(method = lm)+ggtitle("MUT")

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_tcf")]

write.xlsx(merged, file="MUT_significativi_bcat_differenziali_tcf.xlsx")

### Confronto significativi TCF small N2.vs.NE vs significativi BCAT.vs.SCR MUT
tcf_small <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_N2.vs.NE_samples_bcat//MUT_geno_cutoff0.05-N2.vs.NE.deseq2.tsv",
                        quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tcf_small <- tcf_small %>% filter(padj < 0.05)
tcf_small <- tcf_small %>% filter(abs(log2FoldChange) > 0.5849625)
names(tcf_small)[names(tcf_small)=="log2FoldChange"] <- "log2FoldChange_tcf_small_N2.vs.NE"
tcf_small$genes <- rownames(tcf_small)

bcat <- read.table("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_bcat.vs.scr/MUT_geno_cutoff0.05-b2.vs.scr.deseq2.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(bcat)[names(bcat)=="log2FoldChange"] <- "log2FoldChange_bcat.vs.scr"
for (i in rownames(bcat)) {
  if (abs(bcat[i,"log2FoldChange_bcat.vs.scr"]) > 0.5849625 & bcat[i,"padj"] < 0.05) {
    bcat[i,"differenziali_bcat"] <- "YES"
  } else {
    bcat[i,"differenziali_bcat"] <- "NO"
  }
}
bcat$genes <- rownames(bcat)

merged <- merge(bcat, tcf_small, by="genes")

ggplot(merged, aes(y=log2FoldChange_bcat.vs.scr, x=log2FoldChange_tcf_small_N2.vs.NE))+geom_point(aes(color=differenziali_bcat))+geom_smooth(method = lm)+ggtitle("MUT")

merged <- merged[,c("genes","log2FoldChange_bcat.vs.scr", "log2FoldChange_tcf_small_N2.vs.NE", "differenziali_bcat")]

write.xlsx(merged, file="MUT_significativi_tcf_differenziali_bcat.xlsx")
