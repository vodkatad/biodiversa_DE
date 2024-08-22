## sistemare le mutmat con i recuperi

library(tidyverse)
library(openxlsx)

casi_f <- snakemake@input[["casi_df"]]
biob_f <- snakemake@input[["result"]]
san_f <- snakemake@input[["result_san"]]
wes_f <- snakemake@input[["result_wes"]]
fra_f <- snakemake@input[["mut_fra"]]
vbiob_f <- snakemake@input[["vaf_result"]]
vsan_f <- snakemake@input[["vaf_result_san"]]
vwes_f <- snakemake@input[["vaf_result_wes"]]
pbiob_f <- snakemake@input[["prot_result"]]
psan_f <- snakemake@input[["prot_result_san"]]
pwes_f <- snakemake@input[["prot_result_wes"]]
bin_tsv <- snakemake@output[["matrix_bin_mut_tsv"]]
bin_ex <- snakemake@output[["matrix_bin_mut_excel"]]
vaf <- snakemake@output[["mvaf"]]
prot <- snakemake@output[["mprot"]]

#casi_f <- "/scratch/trcanmed/DE_RNASeq/local/share/data/chemio_def_jul23/chemiojul23_Extended_Data_Table.xlsx"
casi <- read.xlsx(casi_f)
casi <- casi$CASE

#biob_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_bin.tsv"
biob <- read.table(biob_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = F)
biob$smodel <- rownames(biob)
biob <- biob %>% filter(smodel %in% casi)
biob$smodel <- NULL
#length(intersect(rownames(biob), casi))

#san_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_bin_sanger.tsv"
san <- read.table(san_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
san$smodel <- rownames(san)
san <- san %>% filter(smodel %in% casi)
san$smodel <- NULL
#length(intersect(rownames(san), casi))

#wes_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_bin_wes.tsv"
wes <- read.table(wes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
wes$smodel <- rownames(wes)
wes <- wes %>% filter(smodel %in% casi)
#length(intersect(rownames(wes), casi))
wes$smodel <- NULL
#wes_model <- rownames(wes)
#san_model <- rownames(san)
#setdiff(wes_model, san_model)

## tutti i dati del wes sono già contenuti nel sanger :( non lo aggiungo a rbind
## vedi pt 186397978

res <- rbind(biob, san)

#duplicated <- res[,duplicated(colnames(res))]
#duplicated <- colnames(duplicated)

# check_casi <- rownames(res)
# length(intersect(casi, check_casi))
# missing <- setdiff(casi, check_casi)
# 
# #fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
# fra <- read.table(fra_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# rownames(fra) <- fra$genes
# fra$genes <- NULL
# tfra <- as.data.frame(t(fra))
# tfra$CASE <- rownames(tfra)
# #rownames(tfra) <- NULL
# 
# tfra <- tfra %>% filter(CASE %in% missing)
# tfra[] <- lapply(tfra, function(x) {
#   x <- gsub("True", 1, x)
#   x <- gsub("False", 0, x)
#   as.integer(x)
# })
# tfra$CASE <- NULL
# 
# setdiff(colnames(res), colnames(tfra))

write.table(res, file=bin_tsv, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.xlsx(res, file=bin_ex)

#vbiob_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_vaf_top10.tsv"
vbiob <- read.table(vbiob_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = F)
vbiob <- as.data.frame(t(vbiob))
vbiob$case <- rownames(vbiob)
vbiob <- vbiob %>% filter(case %in% casi)
vbiob$case <- NULL

#vsan_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_vaf_top10_sanger.tsv"
vsan <- read.table(vsan_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsan <- as.data.frame(t(vsan))

#vwes_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_vaf_top10_wes.tsv"
vwes <- read.table(vwes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vwes <- as.data.frame(t(vwes))

res_vaf <- rbind(vbiob, vsan)
## adding "CRC1472" with all 0
res_vaf <- as.data.frame(t(res_vaf))
res_vaf$CRC1472 <- 0
res_vaf <- as.data.frame(t(res_vaf))

write.xlsx(res_vaf, file=vaf, rowNames = TRUE)
#pbiob_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_protein_top10.tsv"
pbiob <- read.table(pbiob_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = F)
pbiob <- as.data.frame(t(pbiob))
pbiob$case <- rownames(pbiob)
pbiob <- pbiob %>% filter(case %in% casi)
pbiob$case <- NULL

#psan_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_protein_top10_sanger.tsv"
psan <- read.table(psan_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
psan <- as.data.frame(t(psan))

#pwes_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/mutmat_protein_top10_wes.tsv"
pwes <- read.table(pwes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pwes <- as.data.frame(t(pwes))

res_protein <- rbind(pbiob, psan)
## adding "CRC1472" with all WT
res_protein <- as.data.frame(t(res_protein))
res_protein$CRC1472 <- "WT"
res_protein <- as.data.frame(t(res_protein))

write.xlsx(res_protein, file=prot, rowNames = TRUE)