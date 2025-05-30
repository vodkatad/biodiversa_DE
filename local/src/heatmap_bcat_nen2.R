bcat1 <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E",
          "CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2",
          "FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2",
          "HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1",
          "MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4",
          "NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53",
          "WNT1","WNT5B","WNT6", "ASCL2", "LGR5", "NKD1")

bcat <- read_xlsx("/home/mferri/metagene_wnt.xlsx")
bcat <- bcat$Metagene_wnt

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/general/fpkm.tsv.gz"
vsd <- read.table(vsd)
vsd$geni <- rownames(vsd)
vsd <- vsd %>% filter(geni %in% bcat)
#setdiff(bcat, vsd$geni)
vsd$geni <- NULL
vsd <- as.data.frame(t(vsd))
vsd$model <- substr(rownames(vsd), 1, 7)

cocomp <- "/home/mferri/tcf7l2_order_quarti.xlsx"
cocomp <- read.xlsx(cocomp)
cocomp <- cocomp$case

vsd <- vsd %>% filter(model %in% cocomp)
vsd$geni <- substr(rownames(vsd), 9, 10)
vsd$all <- paste0(vsd$model, "_", vsd$geni)
vsd$model <- NULL
vsd$geni <- NULL
vsdnorepliche <- vsd %>% filter(all %in% c("CRC0743_NE", "CRC0743_N2"))
rownames(vsdnorepliche) <- vsdnorepliche$all
vsdnorepliche$all <- NULL
vsdnorepliche <- as.data.frame(t(vsdnorepliche))
vsd <- vsd %>% filter(!all %in% c("CRC0743_NE", "CRC0743_N2"))

unique_pairs <- unique(vsd$all)

# Split the dataframe based on unique pairs
df_list <- lapply(unique_pairs, function(combination) {
  subset_df <-  vsd[grep(combination, vsd$all), ]
  subset_df$all <- NULL
  return(subset_df)
})

t_df_list <- lapply(seq_along(df_list), function(i) {
  transposed_df <- as.data.frame(t(df_list[[i]]))
  return(transposed_df)
})

df_list <- lapply(t_df_list, function(df) {
  nuovo_nome_colonna <- substr(names(df)[1], 1, 10)
  if (ncol(df) == 3) {
    df <- mutate(df, !!nuovo_nome_colonna := rowMeans(df[, 1:3]))
    df <- subset(df, select = -c(1, 2,3))
  } else if (ncol(df) == 2) {
    df <- mutate(df, !!nuovo_nome_colonna := rowMeans(df[, 1:2]))
    df <- subset(df, select = -c(1, 2))
  }
  return(df)
})

nuovo_vsd <- do.call(cbind, df_list)
nuovo_vsd <- cbind(nuovo_vsd, vsdnorepliche)
nuovo_vsd <- log(nuovo_vsd+1)

nuovo_vsd <- as.data.frame(t(nuovo_vsd))
nuovo_vsd$type <- substr(rownames(nuovo_vsd), 9,10)
order <- nuovo_vsd
order$cases <- rownames(order)
#order <- order[,c(46,45)]
order <- order[,c(16,17)]
order_df <- order
order_df <- order_df[order(order_df$type, decreasing=TRUE),]
order_df$cases <- NULL
nuovo_vsd$type <- NULL
nuovo_vsd <- as.data.frame(t(nuovo_vsd))

nuovo_vsd <- nuovo_vsd[,rownames(order_df)]

pheatmap(nuovo_vsd, color = colorRampPalette(c("blue","white","red"))(100),scale='row', cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = order_df)
