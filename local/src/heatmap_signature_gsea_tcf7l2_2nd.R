##Read files named
filenames <- list.files(path="/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/GSEA_C2/",
                        pattern="*.tsv")

##Create list of data frame names without the extra part 
names <-substr(filenames,14,20)

###Load all files
for(i in names){
  filepath <- file.path("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/GSEA_C2/",paste("GSEA_results_", i,"_C2_geno_cutoff0.05-NE.vs.N2.tsv",sep=""))
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

pheatmap(cast_combined, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)

## mettere in un vettore la lista unica di gsea
### lista modelli 
## creare matrice righe e colonne 
## due for sui due vettori e per ogni coppia prendo il NES