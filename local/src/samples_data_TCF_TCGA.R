library(tidyverse)

mut <- snakemake@input[["mut_tcf"]]
wt <- snakemake@input[["wt_tcf"]]
meta <- snakemake@output[["meta"]]

#mut <- "//mnt/trcanmed/snaketree/prj/whatever/dataset/tcf7l2/TCGA_coadread_clinical_tcf7l2-mut.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#wt <- "/mnt/trcanmed/snaketree/prj/whatever/dataset/tcf7l2/TCGA_coadread_clinical_tcf7l2-wt.tsv"
wt <- read.table(wt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mut <- mut[,c("Sample.ID", "TCF7L2")]
mut$geno <- "MUT"
mut$TCF7L2 <- NULL

wt <- wt[,c("Sample.ID", "TCF7L2")]
wt$geno <- "WT"
wt$TCF7L2 <- NULL

ss <- rbind(mut, wt)
ss <- ss %>% filter(!Sample.ID == "TCGA-G4-6317-02")
ss$Sample.ID <- gsub("-01", "", ss$Sample.ID)

ss$Sample.ID <- gsub("-", ".", ss$Sample.ID)

## remove TCGA.AA.3558 perché non c'è nei counts

ss <- ss%>% filter(!Sample.ID == "TCGA.AA.3558")
ss$geno <- as.factor(ss$geno)

write.table(ss, meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

