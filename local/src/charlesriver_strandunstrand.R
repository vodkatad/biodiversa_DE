uns_f <- snakemake@input[['expr_matrix_ustr']]
str_f  <- snakemake@input[['expr_matrix_str']]
meta_f <- snakemake@input[['meta']]
output_f <- snakemake@output[['expr_picked']]


meta <- read.table(meta_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
str <- read.table(gzfile(str_f), sep="\t", header=TRUE, stringsAsFactors=FALSE)
uns <- read.table(gzfile(uns_f), sep="\t", header=TRUE, stringsAsFactors=FALSE)


keep_str <-  meta[meta$stranded == 'firststrand','id']
keep_unstr <-  meta[meta$stranded == 'unstranded','id']
str_expr <- str[, colnames(str) %in% c(keep_str, 'Geneid')]
uns_expr <- str[, colnames(uns) %in% c(keep_unstr, 'Geneid')]

stopifnot(length(intersect(colnames(str_expr),colnames(uns_expr)))==1)

res <- merge(str_expr, uns_expr, by="Geneid")
write.table(res, gzfile(output_f), sep="\t", quote=FALSE, row.names=FALSE)
