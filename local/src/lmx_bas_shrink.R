inputLMX <- snakemake@input[['lmx']]
inputShrink <- snakemake@input[['f1']]
output <- snakemake@output[['f2']]

#read the file with the specific genealogy ids + batch number of the aliquots that were sequenced for the treated models
df_longId <- read.table(inputLMX, sep="\t", header=TRUE)

#subset the data frame keeping the observations of the "LMX_BASALE" type with the id and batch columns only
lmx_bas <- subset(df_longId, type == "LMX_BASALE", c("sample_id_R", "batch"))

#extract the short ids and create a new column to save them
lmx_bas$short_id <- substr(lmx_bas$sample_id_R, 0, 7)

#read the file with the shrinkage data
shrinkage <- read.table(inputShrink, sep="\t", header=TRUE)

#merge lmx_bas and shrinkage data frames together using the short ids
lmx_bas_shrink <- merge(lmx_bas, shrinkage, by.x="short_id", by.y="model")

#keep the "sample_id_R", "batch", and "shrink" columns only
lmx_bas_shrink <- subset(lmx_bas_shrink, select=c("sample_id_R", "batch", "shrink"))

#save the lmx_bas_shrink data frame as a tsv file
write.table(lmx_bas_shrink, file=output, quote=FALSE, sep='\t')
