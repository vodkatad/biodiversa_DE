## recupero casi mancanti che hanno il dato di chemio ma non sono in biobanca
library(tidyverse)

undici_f <- snakemake@input[["odiati11"]]
wes_f <- snakemake@input[["wes_txt"]]
ventisette_f <- snakemake@input[["odiati27"]]
sanger_f <- snakemake@input[["sanger_tsv"]]
res_wes <- snakemake@output[["result_wes"]]
vaf_tsv_wes <- snakemake@output[["vaf_result_wes"]]
prot_tsv_wes <- snakemake@output[["prot_result_wes"]]
res_san <- snakemake@output[["result_san"]]
vaf_tsv_san <- snakemake@output[["vaf_result_san"]]
prot_tsv_san <- snakemake@output[["prot_result_san"]]


## wes velculescu

#undici_f <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/casi_da_recuperare_vWES_11odiati.tsv"
undici <- read.table(undici_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

#wes_f <- "/scratch/trcanmed/pdxopedia/local/share/data/old_seq/DTB_WESdata_CompleteStudy_muts_manualfixes.txt"
wes <- read.table(wes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

geni <- c("TP53", "APC", "MYC", "KRAS", "BRCA2", "GNAS", "ARID1A", "KAT6A", "ERBB2", "PIK3CA", "CTNNB1", "NRAS", "BRAF")

wes11 <- wes %>% filter(CASE %in% undici$V1)
wes11 <- wes11 %>% filter(Gene.Symbol %in% geni)

wes11 <- wes11[,c(2,4,11,16)]
wes11prot <- wes11
wes11prot$X..Mutant.Tags <- NULL
wes11prot <- as.data.frame(wes11prot %>% pivot_wider(names_from = CASE, values_from = Amino.Acid..protein.))
wes11prot <- as.data.frame(lapply(wes11prot, gsub, pattern='c', replacement=''))
wes11prot <- as.data.frame(lapply(wes11prot, gsub, pattern='[()]', replacement=''))
wes11prot <- as.data.frame(lapply(wes11prot, gsub, pattern='"', replacement=''))
wes11prot <- as.data.frame(lapply(wes11prot, gsub, pattern='"', replacement=''))
rownames(wes11prot) <- wes11prot$Gene.Symbol
wes11prot$Gene.Symbol <- NULL
wes11prot <- data.frame(lapply(wes11prot, as.character), stringsAsFactors=FALSE, row.names = rownames(wes11prot))
wes11prot[wes11prot == 'NULL'] <- "WT"
wes11prot[is.na(wes11prot)] <- "NA"

foundgenes <- unique(wes11$Gene.Symbol)
missingwt <- setdiff(geni, foundgenes)

wes11prot <- as.data.frame(t(wes11prot))
wes11prot$MYC <- "WT"
wes11prot$BRCA2 <- "WT"
wes11prot$ARID1A <- "WT"
wes11prot$KAT6A <- "WT"
wes11prot$ERBB2 <- "WT"
wes11prot$CTNNB1 <- "WT"
wes11prot$BRAF <- "WT"
wes11prot <- as.data.frame(t(wes11prot))

write.table(wes11prot, file=prot_tsv_wes, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

wes11vaf <- wes11
wes11vaf$Amino.Acid..protein. <- NULL
wes11vaf <- as.data.frame(wes11vaf %>% pivot_wider(names_from = CASE, values_from = X..Mutant.Tags))
wes11vaf <- as.data.frame(lapply(wes11vaf, gsub, pattern='c', replacement=''))
wes11vaf <- as.data.frame(lapply(wes11vaf, gsub, pattern='[()]', replacement=''))
wes11vaf <- as.data.frame(lapply(wes11vaf, gsub, pattern='"', replacement=''))
wes11vaf <- as.data.frame(lapply(wes11vaf, gsub, pattern='"', replacement=''))
rownames(wes11vaf) <- wes11vaf$Gene.Symbol
wes11vaf$Gene.Symbol <- NULL
wes11vaf <- data.frame(lapply(wes11vaf, as.character), stringsAsFactors=FALSE, row.names = rownames(wes11vaf))
wes11vaf[wes11vaf == 'NULL'] <- 0
wes11vaf[is.na(wes11vaf)] <- "NA"

wes11vaf_2 <- as.data.frame(t(wes11vaf))
wes11vaf_2$MYC <- 0
wes11vaf_2$BRCA2 <- 0
wes11vaf_2$ARID1A <- 0
wes11vaf_2$KAT6A <- 0
wes11vaf_2$ERBB2 <- 0
wes11vaf_2$CTNNB1 <- 0
wes11vaf_2$BRAF <- 0
wes11vaf_2 <- as.data.frame(t(wes11vaf_2))


write.table(wes11vaf_2, file=vaf_tsv_wes, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

wes11binary <- wes11vaf
wes11binary <- data.frame(lapply(wes11binary, as.numeric), stringsAsFactors=FALSE, row.names = rownames(wes11binary))
wes11binary[wes11binary > 0] <- 1 
wes11binary[is.na(wes11binary)] <- 1
wes11binary <- as.data.frame(t(wes11binary))
wes11binary$MYC <- 0
wes11binary$BRCA2 <- 0
wes11binary$ARID1A <- 0
wes11binary$KAT6A <- 0
wes11binary$ERBB2 <- 0
wes11binary$CTNNB1 <- 0
wes11binary$BRAF <- 0

write.table(wes11binary, file =res_wes, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)


## sanger

#ventisette_f <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/casi_da_recuperare_sanger_27odiati.tsv"
ventisette <- read.table(ventisette_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

#sanger_f <- "/scratch/trcanmed/pdxopedia/dataset/sanger_targeted_v2_genealogy/driver_muts_all.tsv"
sanger <- read.table(sanger_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

sanger <- sanger %>% filter(Sample %in% ventisette$V1)
sanger <- sanger %>% filter(Gene %in% geni)
sanger$vaf <- sanger$Subs_VAF
sanger[5, "vaf"] <- 0.992
sanger$vaf <- as.numeric(sanger$vaf)
sanger$CASE <- substr(sanger$Sample, 1, 7)
sanger <- sanger[,c(15, 6,7,14)]

sangerprot <- sanger
sangerprot$vaf <- NULL
sangerprot <- as.data.frame(sangerprot %>% pivot_wider(names_from = CASE, values_from = Protein))
sangerprot <- as.data.frame(lapply(sangerprot, gsub, pattern='c', replacement=''))
sangerprot <- as.data.frame(lapply(sangerprot, gsub, pattern='[()]', replacement=''))
sangerprot <- as.data.frame(lapply(sangerprot, gsub, pattern='"', replacement=''))
sangerprot <- as.data.frame(lapply(sangerprot, gsub, pattern='"', replacement=''))
rownames(sangerprot) <- sangerprot$Gene
sangerprot$Gene <- NULL
sangerprot <- data.frame(lapply(sangerprot, as.character), stringsAsFactors=FALSE, row.names = rownames(sangerprot))
sangerprot[sangerprot == 'NULL'] <- "WT"

foundgenes <- unique(sanger$Gene)
missingwt <- setdiff(geni, foundgenes)

sangerprot <- as.data.frame(t(sangerprot))
sangerprot$MYC <- "WT"
sangerprot$BRCA2 <- "WT"
sangerprot$KAT6A <- "WT"
sangerprot <- as.data.frame(t(sangerprot))

write.table(sangerprot, file=prot_tsv_san, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

sangervaf <- sanger
sangervaf$Protein <- NULL
sangervaf <- as.data.frame(sangervaf %>% pivot_wider(names_from = CASE, values_from = vaf))
sangervaf <- as.data.frame(lapply(sangervaf, gsub, pattern='c', replacement=''))
sangervaf <- as.data.frame(lapply(sangervaf, gsub, pattern='[()]', replacement=''))
sangervaf <- as.data.frame(lapply(sangervaf, gsub, pattern='"', replacement=''))
sangervaf <- as.data.frame(lapply(sangervaf, gsub, pattern='"', replacement=''))
rownames(sangervaf) <- sangervaf$Gene
sangervaf$Gene <- NULL
sangervaf <- data.frame(lapply(sangervaf, as.character), stringsAsFactors=FALSE, row.names = rownames(sangervaf))
sangervaf[sangervaf == 'NULL'] <- 0

sangervaf_2 <- as.data.frame(t(sangervaf))
sangervaf_2$MYC <- 0
sangervaf_2$BRCA2 <- 0
sangervaf_2$KAT6A <- 0
sangervaf_2 <- as.data.frame(t(sangervaf_2))


write.table(sangervaf_2, file=vaf_tsv_san, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

sangerbinary <- sangervaf
sangerbinary <- data.frame(lapply(sangerbinary, as.numeric), stringsAsFactors=FALSE, row.names = rownames(sangerbinary))
sangerbinary[sangerbinary > 0] <- 1 
sangerbinary[is.na(sangerbinary)] <- 1
sangerbinary <- as.data.frame(t(sangerbinary))

sangerbinary$MYC <- 0
sangerbinary$BRCA2 <- 0
sangerbinary$KAT6A <- 0

write.table(sangerbinary, file=res_san, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

print("RICORDA DI MODIFICARE A MANO LE PROTEINE DEL WES E ? DELLE PROTEINE DEL SANGER")
