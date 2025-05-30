matrix_f <- snakemake@input[["matrix"]]
meta_f <- snakemake@input[["meta"]]
hmat_f <- snakemake@output[["hmat"]]

data <- read.table(gzfile(matrix_f), sep="\t", header=TRUE)
rownames(data) <- data$Geneid
data$Geneid <- NULL
colnames(data) <- gsub("X", "", colnames(data))
colnames(data) <- substr(colnames(data), 1, 12)
meta <- read.table(meta_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(meta) <- c("old_id", "new_id")
meta$old_id <- gsub('-','.', meta$old_id)
meta$new_id <- gsub('-','_', meta$new_id)
meta <- meta[match(colnames(data),meta$old_id),]
if (!all(colnames(data) == meta$old_id)) {
  stop('Cosa fai, llama?')
}

colnames(data) <- meta$new_id
## elimino momentaneamnte CRC0152LMO0C01003001VT0700R e CRC0152LMO0C04010001R01000 in attesa della nuova analisi
#data$CRC0152LMO0C04010001R01000 <- NULL
#data$CRC0152LMO0C01003001VT0700R <- NULL

#elimino il CRC1961LMO0A01003001VT0500R e CRC1961LMO0A02008001VT0300R per failed strandness
#data$CRC1961LMO0A01003001VT0500R<- NULL
#data$CRC1961LMO0A02008001VT0300R <- NULL

write.table(data, gzfile(hmat_f), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)