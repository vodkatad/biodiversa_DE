
expr_input <- snakemake@input[["expr"]]
metadata_input <- snakemake@input[["metadata"]]
buoni <- snakemake@input[["goodOrBad"]]
expr_output <- snakemake@output[["expr"]]
metadata_output <- snakemake@output[["metadata_clu"]]

save.image('pippo.Rdata')
selection <- read.table(buoni, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote='')
metadata <- read.table(metadata_input, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote='')

#expr_df <- read.table(gzfile(expr_input), sep="\t", header=TRUE, stringsAsFactors=FALSE)
expr_df <- read.table(expr_input, sep="\t", header=TRUE, stringsAsFactors=FALSE)

good <- selection[selection$buoni, 'smodel']

lmo_good <- metadata[metadata$sample %in% good & metadata$type=="LMO_BASALE", ]

lmx_good <- metadata[metadata$sample %in% lmo_good$sample & metadata$type=="LMX_BASALE", ]
lmh_good <- metadata[metadata$sample %in% lmo_good$sample & metadata$type=="LMH", ]

has_h_has_x <- c(lmx_good$sample, lmh_good$sample)

lmo_good_has_xh <- lmo_good[lmo_good$sample %in% has_h_has_x,]

wanted <- rbind(lmo_good_has_xh, lmx_good, lmh_good)

expr_selected <- expr_df[,colnames(expr_df) %in% rownames(wanted)]

write.table(expr_selected, gzfile(expr_output), sep="\t", quote=FALSE)
write.table(wanted, metadata_output, sep="\t", quote=FALSE)


