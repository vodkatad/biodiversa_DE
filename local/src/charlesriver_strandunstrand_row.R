uns_f <- snakemake@input[['matrix_ustr']]
str_f  <- snakemake@input[['matrix_str']]
meta_f <- snakemake@input[['meta']]
output_f <- snakemake@output[['picked']]


meta <- read.table(meta_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
str <- read.table(str_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
uns <- read.table(uns_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)


keep_str <-  meta[meta$stranded == 'firststrand','id']
keep_unstr <-  meta[meta$stranded == 'unstranded','id']
str_reads <- str[str$sample %in% keep_str,]
uns_reads <- uns[uns$sample %in% keep_unstr,]

stopifnot(length(intersect(str_reads$sample,uns_reads$sample))==0)

res <- rbind(str_reads, uns_reads)
write.table(res, output_f, sep="\t", quote=FALSE, row.names=FALSE)
