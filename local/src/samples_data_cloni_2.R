library(stringr)


gep <- snakemake@input[["samples_data"]]

res <- snakemake@output[["tsv"]]

data <- read.table(gzfile(gep), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(data)
data<-data[!grepl("CL",data$geno),]
data<-data[!grepl("3",data$replica),]
data<-data[grep("Cas9",data$geno),]
write.table(data, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)