## bcat

bcat <- "/home/mferri/GO0016055_wnt_signaling_pathway.tsv"
bcat <- read.table(bcat, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bcat <- unique(bcat$SYMBOL)

deg <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/geno_cutoff0.05-b2.vs.scr.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)

deg <- deg %>% filter(genes %in% bcat)

d <- deg[,c("genes", "log2FoldChange")]
d <- d[order(d$log2FoldChange, decreasing = TRUE),]
d$genes <- NULL
d <- as.matrix(d)

minv <- -4
maxv <- 4
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
         breaks = bk, color=my_palette,  na_col = "#FFFFFF")

##jak

jak <- "/home/mferri/GO0007259_pathway_JAK-STAT.tsv"
jak <- read.table(jak, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
jak <- unique(jak$SYMBOL)

deg <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/bcat.vs.scr/geno_cutoff0.05-b2.vs.scr.deseq2.tsv"
deg <- read.table(deg, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
deg$genes <- rownames(deg)

deg <- deg %>% filter(genes %in% jak)

d <- deg[,c("genes", "log2FoldChange")]
d <- d[order(d$log2FoldChange, decreasing = TRUE),]
d$genes <- NULL
d <- as.matrix(d)

minv <- -2
maxv <- 2
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
         breaks = bk, color=my_palette,  na_col = "#FFFFFF")



