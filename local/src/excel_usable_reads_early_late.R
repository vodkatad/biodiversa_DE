## usable reads biobanca early late

library(tidyverse)
library(stringr)
library(openxlsx)

df_f <- snakemake@input[["reads"]]
nomi_f <- snakemake@input[["meta"]]
id_f <- snakemake@input[["ega"]]
alias_f <- snakemake@input[["run"]]
excel <- snakemake@output[["us_reads"]]

#df_f <- "/mnt/cold1/bioinfotree/prj/ngs_standard_processing/dataset/earlylate_biobanca/usable_reads"
df <- read.table(df_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df <- df[, c(1,2,7)]
df$fraction <- df$Assigned/df$tot
df$sample <- sapply(str_split(df$sample, "_",  n = 3), `[`, 1)
df$sample <- gsub("V1", "", df$sample)

#nomi_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/early_late/ss_early_late_biobanca.txt"
nomi <- read.table(nomi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(nomi) <- c("sample", "nomi")

df <- merge(df, nomi, by="sample")
df$sample <- NULL
df <- df[,c(4,1,2,3)]
colnames(df) <- c("Genealogy ID",	"Total reads",	"Total mapped reads",	"Fraction of mapped reads")

#id_f <- "/mnt/cold2/snaketree/prj/EGAs/earlylate/EGAsubmitter/dataset/submission/EGASTORE/samples_EGAID.tsv"
id <- read.table(id_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(id) <- c("Genealogy ID", "EGA Sample Accession ID")

#setdiff(id$`Genealogy ID`, df$`Genealogy ID`)
#[1] "CRC0152LMO0C04010001R01000-1"

df[41, "Genealogy ID"] <- "CRC0152LMO0C04010001R01000-1"

df <- merge(df, id, by="Genealogy ID")

#alias_f <- "/mnt/cold2/snaketree/prj/EGAs/earlylate/EGAsubmitter/dataset/submission/EGASTORE/runs_EGAID.tsv"
alias <- read.table(alias_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
alias <- alias[,c(2,7)]
alias$sample <- sapply(str_split(alias$sample, ",",  n = 3), `[`, 2)
alias$sample <- gsub(" 'alias': ", "", alias$sample)
alias$sample <- gsub("'", "", alias$sample)
colnames(alias) <- c("EGA Run Accession ID", "Genealogy ID")

df <- merge(df, alias, by="Genealogy ID")

write.xlsx(df, excel)