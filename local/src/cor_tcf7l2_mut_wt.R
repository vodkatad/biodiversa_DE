## correlation tcf7 

mut <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/MUT_NE.vs.N2/MUT_geno_cutoff0.05-NE.vs.N2.deseq2.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

wt <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/WT_geno_cutoff0.05-NE.vs.N2.deseq2.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# means_sds <- function(df) {
#   means <- apply(df, 1, mean)
#   med <- median(means)
#   df <- as.data.frame(df[means > med,])
#   sds <- apply(df, 1, sd) 
#   sds <- sds[order(-sds)]
#   n <- length(sds)
#   keep <- head(sds, round(0.10*n))
#   keep_genes <- names(keep)
#   df <- df[rownames(df) %in% keep_genes,]
#   df$genes <- rownames(df)
#   df<- df[,c(7,2)]
#   return(df)
# }
# 
# mut <- means_sds(mut)
# wt <- means_sds(wt)
mut$log2FoldChange <- mut$log2FoldChange*-1
wt$log2FoldChange <- wt$log2FoldChange*-1

names(mut)[names(mut)=="log2FoldChange"] <- "MUT_LFC"
names(wt)[names(wt)=="log2FoldChange"] <- "WT_LFC"
names(mut)[names(mut)=="padj"] <- "padj_MUT"
names(wt)[names(wt)=="padj"] <- "padj_WT"

mut$genes <- rownames(mut)
wt$genes <- rownames(wt)

merged <- merge(mut, wt, by="genes")
merged <- merged[,c("genes", "MUT_LFC", "WT_LFC", "padj_MUT", "padj_WT")]

bcat_m <- merged

int <- c("ADAR","AZI2","CACTIN","CDC37","CH25H","CNOT7","DCST1","EIF4E2",
         "FADD","GIGYF2","HDAC4","IFI27","IFIH1","IFIT1","IFITM1","IFITM2","IFITM3","IFNA1","IFNA10",
         "IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8",
         "IFNAR1","IFNAR2","IFNB1","IFNE","IFNK","IFNW1","IKBKE","IRAK1","IRF3","IRF7","ISG15","JAK1",
         "LSM14A","MAVS","METTL3","MIR21","MMP12","MUL1","MX1","MYD88","NLRC5","OAS1","OAS2","OAS3","PTPN1",
         "PTPN11","PTPN2","PTPN6","RBM47","RNF185","SAMHD1","SETD2","SHFL","SHMT2","SIN3A","SMPD1","SP100",
         "STAT1","STAT2","STING1","TANK","TBK1","TBKBP1","TRAF3","TREX1","TRIM41","TRIM56","TRIM6","TRIM65",
         "TTLL12","TYK2","UBE2K","USP18","USP27X","USP29","WNT5A","YTHDF2","YTHDF3","ZBP1")

rownames(merged) <- merged$genes
for (i in rownames(merged)) {
  if (i %in% int) { 
    merged[i, "int"] <- "yes"
  } else {
    merged[i, "int"] <- "no"
  }
}

merged$color <- ifelse(merged$int == "yes", "yes_color", "no_color")

yes_color <- "red"
no_color <- "black"

ggplot(merged, aes(x=MUT_LFC, y=WT_LFC, color=color))+geom_point()+
  scale_color_manual(values = c(yes_color = yes_color, no_color = no_color))

int <- merged %>% filter(int == "yes")

ggplot(int, aes(x=MUT_LFC, y=WT_LFC))+geom_point()+geom_smooth(method = lm)
cor.test(int$MUT_LFC, int$WT_LFC)

cor.test(merged$MUT_LFC, merged$WT_LFC)

# bcat <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E",
#            "CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2",
#            "FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2",
#            "HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1",
#            "MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4",
#            "NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53",
#            "WNT1","WNT5B","WNT6", "ASCL2", "LGR5", "NKD1")

# bcat <- read.xlsx("/home/mferri/Diret_targets_Wnt_homepage.xlsx")
# bcat <- bcat$Gene

bcat <- read.xlsx("/home/mferri/FINAL_WNT_TARGETS.xlsx", colNames = FALSE)
bcat <- bcat$X1

bcat_m$padj_MUT <- NULL
bcat_m$padj_WT <- NULL

bcat_m <- bcat_m %>% filter(genes %in% bcat)
rownames(bcat_m) <- bcat_m$genes
bcat_m$genes <- NULL

# row_means <- rowMeans(bcat_m)
# ordered_row_names <- rownames(bcat_m)[order(row_means, decreasing = FALSE)]
# bcat_m <- bcat_m[ordered_row_names, ]
# bcat_m$padj_MUT <- NULL
# bcat_m$padj_WT <- NULL

ordered_bcat <- bcat_m[order(bcat_m$WT_LFC),]
#ordered_bcat_delta <- bcat_m
#ordered_bcat_delta$delta <- sign(ordered_bcat_delta$WT_LFC)*abs((abs(ordered_bcat_delta$WT_LFC) - abs(ordered_bcat_delta$MUT_LFC)))
#ordered_bcat_delta <- ordered_bcat_delta[order(ordered_bcat_delta$delta),]

pheatmap(bcat_m, cluster_rows = FALSE, cluster_cols = FALSE)

d <- ordered_bcat

minv <- -1
maxv <- 1
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#e1e1e1", "#e1e1e1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette)

