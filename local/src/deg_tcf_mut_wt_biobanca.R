library(tidyverse)

meda <- snakemake@input[["metadati_biobanca"]]
new_f <- snakemake@input[["new_seq"]]
pdx_f <- snakemake@input[["appello"]] 
meta <- snakemake@output[["meta"]]

tcfmut <- c("CRC0059", "CRC0065", "CRC0068", "CRC0076", "CRC0123",
      "CRC0148","CRC0152", "CRC0169", "CRC0173", "CRC0196",
      "CRC0204", "CRC0277", "CRC0316", "CRC0324", "CRC0327",
      "CRC0355", "CRC0399", "CRC0416", "CRC0427", "CRC0449",
      "CRC0464", "CRC0516", "CRC0517", "CRC0556", "CRC0556",
      "CRC0591", "CRC1063", "CRC1067", "CRC1138", "CRC1182",
      "CRC1239", "CRC1245", "CRC1278", "CRC1306", "CRC1331",
      "CRC1336", "CRC1342", "CRC1432", "CRC1473", "CRC1474",
      "CRC1477", "CRC1586", "CRC1598", "CRC1629", "CRC1629",
      "CRC1675", "CRC1709", "CRC1729", "CRC1729", "CRC1811",
      "CRC1917", "CRC1963", "CRC2236", "CRC2247", "CRC2248",
      "CRC2307", "CRC2307", "CRC2378", "CRC2379", "CRC2384",
      "CRC2384", "CRC2385", "CRC2388", "CRC2403", "CRC2403",
      "CRC2546", "CRC2734", "CRC2870", "CRC3196", "CRC3196")

#meda <- "/mnt/cold1//snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"

meda <- read.table(meda, sep='\t', quote="", header=TRUE)
meda$sample_id <- gsub('-', '.', meda$sample_id_R, fixed = TRUE)
meda <- meda %>%
  filter(type %in% c("LMX_BASALE", "LMX_BASALE.1"))
meda$ssample <- substr(meda$sample_id, 1, 7)
tcfmutrna<- meda %>% filter(ssample %in% tcfmut)

#new_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/chemio_collection/samples_data"
new <- read.table(new_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
new <- new %>% filter(batch == "new")  
new <- new %>% filter(id %in% tcfmut)  

# xeno in qualsiasi appello seq e ha un basale
# escludendo i bad boys 
# prendo 40 wt a caso che hanno il basale

#pdx_f <- "/home/mferri/appello_sequenziamenti_PDX.tsv"
pdx <- read.table(pdx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)  
pdx <- pdx %>% filter(anyseq == "yes")
pdx <- pdx %>% filter(!smodel %in% tcfmut)
tcfwt <- pdx$smodel

tcfwtrna <- meda %>% filter(ssample %in% tcfwt)
set.seed(123)
random_strings <- sample(unique(tcfwtrna$ssample), 37)

tcfwtrna <- tcfwtrna %>% filter(ssample %in% random_strings)
tcfwtrna$geno <- "WT"
tcfmutrna$geno <- "MUT"

tcfwtrna <- tcfwtrna[,c("sample_id", "ssample", "batch", "geno")]
tcfmutrna <- tcfmutrna[,c("sample_id", "ssample", "batch", "geno")]

tcfbiob <- rbind(tcfwtrna, tcfmutrna)
names(tcfbiob)[names(tcfbiob)=="ssample"] <- "model"
rownames(tcfbiob) <- tcfbiob$sample_id
tcfbiob$sample_id <- NULL

write.table(tcfbiob, file=meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)