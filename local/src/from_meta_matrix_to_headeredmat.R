matrix_f <- snakemake@input[["matrix"]]
meta_f <- snakemake@input[["meta"]]
hmat_f <- snakemake@output[["hmat"]]

data <- read.table(gzfile(matrix_f), sep="\t", header=TRUE, row.names=1)
meta <- read.table(meta_f, sep="\t", header=TRUE)

meta <- meta[match(colnames(data),meta$sample_id_R),]
if (!all(colnames(data) == meta$sample_id_R)) {
    stop('Cosa fai, llama?')
}

meta$new_id <- paste0(substr(meta$sample_id_R,0,7),'_',meta$type)
colnames(data) <- meta$new_id

write.table(data, gzfile(hmat_f), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)