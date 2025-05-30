##Read files named
filenames <- list.files(path="/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/GSEA_H/",
                        pattern="*.tsv")

##Create list of data frame names without the extra part 
names <-substr(filenames,14,20)

###Load all files
for(i in names){
  filepath <- file.path("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/TCF7L2_DEG/GSEA_H/",paste("GSEA_results_", i,"_H_geno_cutoff0.05-NE.vs.N2.tsv",sep=""))
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
#[1] "CRC0743"
order_dep$CRC0743 <- NULL
#order_dep$CRC0080 <- NULL
order_dep <- as.data.frame(t(order_dep))
order_dep <- order_dep[order(order_dep$mut),]

cast_combined <- as.data.frame(cast_combined)
cast_combined_ordered <- as.data.frame(t(cast_combined))
cast_combined_ordered <- cast_combined[, rownames(order_dep)]


pheatmap(cast_combined_ordered, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, annotation_col = order_dep, fontsize_row = 5)

