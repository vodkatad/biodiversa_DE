matrix_f <- snakemake@input[["matrix"]]
meta_f <- snakemake@input[["meta"]]
hmat_f <- snakemake@output[["hmat"]]

data <- read.table(gzfile(matrix_f), sep="\t", header=TRUE, row.names=1)
meta <- read.table(meta_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
meta$old_id <- gsub('-','.', meta$old_id)
meta$old_id <- paste0('X', meta$old_id)
colnames(data) <- sapply(strsplit(colnames(data), '_'), function(x) {x[[1]][1]})
meta <- meta[match(colnames(data),meta$old_id),]
if (!all(colnames(data) == meta$old_id)) {
    stop('Cosa fai, llama?')
}

colnames(data) <- meta$new_id

write.table(data, gzfile(hmat_f), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
