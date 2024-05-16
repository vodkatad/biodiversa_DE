counts_1f <- snakemake@input[["counts"]]
counts_2f <- snakemake@input[["add"]]
meta_f <- snakemake@output[["meta"]]
ocounts_f <- snakemake@output[["counts"]]

counts_1 <- read.table(gzfile(counts_1f), sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
counts_2 <- read.table(gzfile(counts_2f), sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)

counts_2 <- counts_2[,!grepl('CAS9', colnames(counts_2))]
counts_2 <- counts_2[,!grepl('CL', colnames(counts_2))]
m <- merge(counts_1, counts_2, by="row.names")
colnames(m) <- gsub('NO_EGF', 'NOEGF', colnames(m))
rownames(m) <- m$Row.names
m$Row.names <- NULL

samples <- data.frame(id=colnames(m), stringsAsFactors=FALSE)
samples$smodel <- substr(samples$id,0,7)
samples$id <- gsub('NO_EGF', 'NOEGF', samples$id)
samples$sorted <- ifelse(grepl('NEG', samples$id), 'NEG', ifelse(grepl('POS', samples$id), 'POS', 'BULK'))
samples$treat <- sapply(strsplit(samples$id, "_"), function(x) {x[[3]]} )
write.table(m, gzfile(ocounts_f), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(samples, meta_f, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
