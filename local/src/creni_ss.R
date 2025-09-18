library(stringr)


gep <- snakemake@input[["meta"]]

res <- snakemake@output[["tsv"]]

data <- read.table(gzfile(gep), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
combo<-data[grepl("Combo_72h",data$trattamento),]
cet<-data[grepl("Cetux_72h",data$trattamento),]
data <- rbind(combo, cet)
write.table(data, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)