vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/b2_general/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- as.data.frame(t(vsd))
vsd$geno <- substr(rownames(vsd), 12, 13)
vsd <- vsd %>% filter(!geno == "t2")
vsd$geno <- NULL
rownames(vsd) <- gsub("sh_", "", rownames(vsd))
vsd <- as.data.frame(t(vsd))
vsd$gene <- rownames(vsd)
bcat <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E",
          "CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2",
          "FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2",
          "HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1",
          "MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4",
          "NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53",
          "WNT1","WNT5B","WNT6", "ASCL2", "LGR5", "NKD1")
#bcat <- read_xlsx("/home/mferri/metagene_wnt.xlsx")
#bcat <- bcat$Metagene_wnt
vsd <- vsd %>% filter(gene %in% bcat)
vsd$gene <- NULL

vsd <- as.data.frame(t(vsd))
vsd$all <- substr(rownames(vsd), 1, 10)

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
nuovo_vsd <- log(nuovo_vsd+1)

order_dep_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx"
order_dep <- read_xlsx(order_dep_f)
order_dep$quartile <- ntile(order_dep$co_comp_N2, 4)
names(order_dep)[names(order_dep)=="case"] <- "model"
order_dep <- as.data.frame(order_dep)
rownames(order_dep) <- order_dep$model

mut <- c("CRC0148", "CRC1331", "CRC0399", "CRC1278", "CRC0277", "CRC0327",
         "CRC1729", "CRC0152", "CRC0196", "CRC0059", "CRC0065", "CRC0464",
         "CRC0316", "CRC1239")
for (i in rownames(order_dep)) {
  if (order_dep[i, "model"] %in% mut) {
    order_dep[i, "mut"] <- "MUT"
  } else {
    order_dep[i, "mut"] <- "WT"
  }
}

order_dep$co_comp_N2 <- NULL
order_dep$model <- NULL
casi <- setdiff(rownames(order_dep), substr(colnames(nuovo_vsd), 1, 7))
order_dep$cases <- rownames(order_dep)
order_dep <- order_dep %>% filter(!cases %in% casi)
order_dep <- rbind(order_dep, order_dep)
order_dep[1:6, "geno"] <- "_sc"
order_dep[7:12, "geno"] <- "_b2"
order_dep$new <- paste0(order_dep$cases, order_dep$geno)
rownames(order_dep) <- order_dep$new
order_dep <- order_dep[,c(1:2)]
order_dep$type <- substr(rownames(order_dep), 9, 10)
order_dep$quartile <- as.factor(order_dep$quartile)

nuovo_vsd_or <- nuovo_vsd[,rownames(order_dep)]

pheatmap(nuovo_vsd, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = order_dep)
pheatmap(nuovo_vsd_or,color = colorRampPalette(c("blue","white","red"))(100),scale='row',cluster_rows = FALSE,cluster_cols = FALSE, annotation_col = order_dep)
