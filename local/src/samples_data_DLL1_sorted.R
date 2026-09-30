output_dir <- snakemake@output[['tsv']]
input_dir <- snakemake@input[['data']]


check_count <- read.table(
  gzfile(input_dir),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
)

samples_keep <- grep("CRC0322",#|^CRC0327_APC
                     colnames(check_count),
                     value = TRUE)


count_subset <- check_count[, samples_keep]

samples_data <- data.frame(
  row.names = samples_keep,
  sample = sub("_.*", "", samples_keep),
  group = sub(".*APC", "", samples_keep),
  stringsAsFactors = FALSE
)

write.table(
  samples_data,
  file = output_dir,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)