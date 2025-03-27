library(stringr)


gep <- snakemake@input[["meta"]]
#ss<-read.table(gep,header=TRUE,row.names = 1,stringsAsFactors = FALSE)
res <- snakemake@output[["tsv"]]

data <- read.table(gzfile(gep), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames_list <- colnames(data)[-1]
split_names <- strsplit(colnames_list, "_")

# Creazione del nuovo dataframe
metadata <- data.frame(
  id = colnames_list,
  sample = sapply(split_names, function(x) x[1]),  # "CRC0322"
  geno = sapply(split_names, function(x) x[2]),    # "CL23", "CL88", "Cas9", "wt"
  trattamento = sapply(split_names, function(x) x[3]), # "CETUX", "EGF0.1"
  replica = sapply(split_names, function(x) gsub("\\D", "", x[4])) # replicati
)
rownames(metadata)<-metadata$id
metadata$id<-NULL

write.table(metadata, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
