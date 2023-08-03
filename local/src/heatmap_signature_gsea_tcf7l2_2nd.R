##Read files named
filenames <- list.files(path="/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/GSEA_C2/",
                        pattern="*.tsv")

##Create list of data frame names without the extra part 
names <-substr(filenames,14,20)

###Load all files
for(i in names){
  filepath <- file.path("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/GSEA_C2/",paste("GSEA_results_", i,"_C2_geno_cutoff0.05-NE.vs.N2.tsv",sep=""))
  assign(i, read.table(filepath, quote="", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
}  

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

add_name_column <- function(df, name) {
  df$Name <- name
  return(df)
}

# Get the names of the data frames in the list
df_names <- names(dfs)

# Iterate over the list and add the name column to each data frame
df_list_with_name <- Map(add_name_column, dfs, df_names)


get_sign <- function(df){
  df <- df %>% filter(p.adjust < 0.05)
}

result_list <- lapply(df_list_with_name, get_sign)
combined_df <- bind_rows(result_list)

cast_combined <- cast(combined_df, ID~Name, value = "NES") 
cast_combined[is.na(cast_combined)] <- 0 
rownames(cast_combined) <- cast_combined$ID
cast_combined$ID <- NULL

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
order_dep <- as.data.frame(t(order_dep))
setdiff(colnames(order_dep), colnames(cast_combined))
#[1] "CRC0743" "CRC0080"
order_dep$CRC0743 <- NULL
order_dep$CRC0080 <- NULL
order_dep <- as.data.frame(t(order_dep))

pheatmap(cast_combined, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, annotation_col = order_dep)

minv <- min(cast_combined)
maxv <- max(cast_combined)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pheatmap(cast_combined, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, annotation_col = order_dep, breaks = bk, color = my_palette)

analysis <- cast_combined

quart12 <- order_dep %>% filter(quartile == 1 | quartile == 2)
quart12 <- rownames(quart12)
analysis <- analysis[, quart12]

analysis$mean <- rowMeans(analysis)
analysis <- analysis %>% filter(mean > 2.5)
analysis$mean <- NULL

selected <- rownames(analysis)
selected_sign <- cast_combined
selected_sign$signature <- rownames(selected_sign)
selected_sign <- selected_sign %>% filter(signature %in% selected)
selected_sign$signature <- NULL
selected_sign <- as.data.frame(selected_sign)

pheatmap(selected_sign, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, annotation_col = order_dep)

write.table(analysis, file = "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/signature_quartili1.2.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
