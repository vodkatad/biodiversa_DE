## heatmap chemio real time

df <- "/home/mferri/CHEMIO_SIMO_R_VS_S_ANALISI_MF.xlsx"
df <- as.data.frame(read_xlsx(df, sheet = "MF"))
df <- df[, !names(df) %in%  c("FBXO18", "RING1", "RAD51", "POLD1")]
df$N <- NULL
rownames(df) <- df$CRC
annot <- df[,c(1,2)]
annot <- annot[order(annot$type, decreasing=TRUE),]
annot$CRC <- NULL
df$CRC <- NULL

order <- df
order <- order %>% filter(type == "Sensitive")
order$type <- NULL
order <- as.data.frame(t(order))
order$mean <- rowMeans(order)
order <- order[order(order$mean, decreasing=TRUE),]

df$type <- NULL

df <- df[rownames(annot), ]
df <- df[,rownames(order)]

pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = annot)

ann_colors = list(
  type = c(Resistant=rgb(red=165,green=0,blue=25, max = 255), Sensitive=rgb(30, 85, 130, max = 255)))
annot$type <- as.factor(annot$type)

minv <- min(df)
maxv <- max(df)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = annot, annotation_colors = ann_colors,breaks = bk, color = my_palette)

minv <- -15
maxv <- 15
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = annot, annotation_colors = ann_colors,breaks = bk, color = my_palette)


