### samples data for deg early late

library(tidyverse)
library(readxl)

df_f <- snakemake@input[["meta"]]
res <- snakemake@output[["tsv"]]

## report
#df_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/early_late/Sequencing_report_biobanca_earlylate.xlsx"
df <- read_xlsx(df_f)
df <- df[c(2,1)]
df$Sample <- NULL
df$Customer_Sample_ID <- gsub('-','_', df$Customer_Sample_ID)
colnames(df) <- c("id")
df$model <- substr(df$id, 1, 7)
df$passage <- substr(df$id, 16, 17)
df$passage <- as.numeric(df$passage)
for (i in seq(rownames(df))) {
  if (df[i,"passage"] == 3) {
    df[i, "passage_el"] <- "early"
  } else {
    df[i, "passage_el"] <- "late"
  }
}
df$passage <- df$passage_el
df$passage_el <- NULL
## elimino momentaneamnte CRC0152LMO0C01003001VT0700R e CRC0152LMO0C04010001R01000 in attesa della nuova analisi
#df <- df[-7,]
#df <- df[-40,]
#elimino il CRC1961LMO0A01003001VT0500R e CRC1961LMO0A02008001VT0300R per failed strandness
#df <- df[-22,]
#df <- df[-22,]
write.table(df, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

