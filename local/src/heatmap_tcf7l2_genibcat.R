vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/general/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- as.data.frame(t(vsd))
vsd$geno <- substr(rownames(vsd), 9, 10)
vsd <- vsd %>% filter(!geno == "N2")
vsd$geno <- NULL

# bcat <- read.xlsx("/home/mferri/supp_gkq363.xlsx", colNames = FALSE)
# bcat <- bcat$X1
# bcat <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E",
#           "CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2",
#           "FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2",
#           "HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1",
#           "MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4",
#           "NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53",
#           "WNT1","WNT5B","WNT6", "ASCL2", "LGR5", "NKD1")

bcat <- read.xlsx("/home/mferri/metagene_wnt.xlsx")
bcat <- bcat$Metagene_wnt

vsd <- as.data.frame(t(vsd))
vsd$gene <- rownames(vsd)
vsd <- vsd %>% filter(gene %in% bcat)
vsd$gene <- NULL
vsd <- as.data.frame(t(vsd))
vsd$all <- substr(rownames(vsd), 1, 7)

vsdnorepliche <- vsd %>% filter(all %in% c("CRC0743"))
rownames(vsdnorepliche) <- vsdnorepliche$all
vsdnorepliche$all <- NULL
vsdnorepliche <- as.data.frame(t(vsdnorepliche))
vsd <- vsd %>% filter(!all %in% c("CRC0743"))
colnames(vsdnorepliche) <- c("CRC0743_NE")

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

setdiff(bcat, rownames(nuovo_vsd))
## all the diff are only genes with another name in the excel

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
#order_dep <- as.data.frame(t(order_dep))
#setdiff(colnames(order_dep), substr(colnames(nuovo_vsd), 1, 7))
#order_dep <- order_dep[, c("CRC0148", "CRC0291", "CRC0322", "CRC0542", "CRC1278", "CRC1430")]
#[1] "CRC0743"
#order_dep$CRC0743 <- NULL
#order_dep$CRC0080 <- NULL
#order_dep <- as.data.frame(t(order_dep))
order_dep <- order_dep[order(order_dep$mut),]
rownames(order_dep) <- paste0(rownames(order_dep), "_NE")
#order_dep <- as.data.frame(t(order_dep))
order_dep <- order_dep[order(order_dep$mut),]
#order_dep <- as.data.frame(t(order_dep))
order_dep$quartile <- as.factor(order_dep$quartile)

nuovo_vsd_ordered <- nuovo_vsd[,rownames(order_dep)]

pheatmap(nuovo_vsd_ordered, color = colorRampPalette(c("blue","white","red"))(100),scale='row', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 5, annotation_col = order_dep)
