matrix_f <- snakemake@input[["matrix"]]
meta_f <- snakemake@input[["meta"]]
hmat_f <- snakemake@output[["hmat"]]
meta_o_f <- snakemake@output[["meta"]]

data <- read.table(gzfile(matrix_f), sep="\t", header=TRUE, row.names=1)
meta <- read.table(meta_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(meta) <- c('id', 'sample')
all_samples <- unique(meta$sample)
meta$id <- paste0('X', meta$id)

sum_samples_lanes <- function(sample, meta, data) {
    wanted <- meta[meta$sample==sample, 'id']
    res <- apply(data[,colnames(data) %in% wanted], 1, sum)
    return(res)
}

sum_data <- sapply(all_samples, sum_samples_lanes, meta, data)

res_meta <- data.frame(id=colnames(sum_data), stringsAsFactors=FALSE)
res_meta$treat <- sapply(strsplit(res_meta$id, "_", fixed=TRUE), function(x) {x[3]})
res_meta$smodel <- sapply(strsplit(res_meta$id, "_", fixed=TRUE), function(x) {x[1]})
res_meta[is.na(res_meta$treat), 'treat'] <- 'EGF'

write.table(sum_data, gzfile(hmat_f), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(res_meta, meta_o_f, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
