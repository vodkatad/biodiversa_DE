### st4 expr dei geni HR nei casi sottoposti ad analisi mutazionale
library(tidyverse)

geni_f <- snakemake@input[["geni_hr"]]
appello_f <- snakemake@input[["casi"]]
vsd_f <- snakemake@input[["expr"]]
res <- snakemake@output[["result"]]
log_f <- snakemake@log[['log']]

#geni_f <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/targeted_homologous.tsv"
genes <- read.table(geni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- genes$gene_symbol

#appello_f <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/appello_mut.tsv"
appello <- read.table(appello_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
appello <- appello$model

#vsd_f <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(vsd) <- gsub("H_", "", rownames(vsd))
vsd_geni <- vsd
vsd_geni$geni <- rownames(vsd_geni)
vsd_geni <- vsd_geni %>% filter(geni %in% genes)
vsd_geni$geni <- NULL
vsd_gc <- as.data.frame(t(vsd_geni))

vsd_gc$smodel <- substr(rownames(vsd_gc), 1, 7)
vsd_gc$type <- substr(rownames(vsd_gc), 8, 10)
vsd_gc <- vsd_gc %>% filter(type == "LMX")
vsd_gc$type <- NULL

vsd_gc <- vsd_gc %>% filter(smodel %in% appello)
casi_appellati <- unique(vsd_gc$smodel)
sink(log_f, append = TRUE)
print("Casi con RNA mancante in biobanca")
setdiff(appello, casi_appellati)
sink()
vsd_gc$smodel <- NULL

write.table(vsd_gc, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
