### preparation clinical data chemio jul23 for circos

library(openxlsx)
library(tidyverse)

#order for circos mail livio 31/08
# Sidedness
# Stage
# Therapy before
# MSI/MSS da recuperare da quello per biobanca
# KRAS da recuperare da quello per biobanca
# NRAS da recuperare da quello per biobanca
# BRAF da recuperare da quello per biobanca
# Sex
# Age at collection

data_f <- "/scratch/trcanmed/DE_RNASeq/local/share/data/chemio_def_jul23/chemiojul23_Extended_Data_Table.xlsx"
data <- read.xlsx(data_f)
data <- as.data.frame(cbind(data$CASE, data$SITE.OF.PRIMARY, data$STAGE, data$`THERAPY.BEFORE.COLLECTION.(Y/N)`, data$SEX, data$`AGE.AT.COLLECTION.(years)`))
colnames(data) <- c("Model", "Sidedness", "Stage", "Therapy_Before", "Sex", "Age_at_Collection")

fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
fra <- read.table(fra_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(fra) <- fra$genes
fra$genes <- NULL
tfra <- as.data.frame(t(fra))
tfra$Model <- rownames(tfra)
rownames(tfra) <- NULL

merged <- merge(data, tfra, by = "Model", all.x = TRUE)

msis_f <- "/scratch/trcanmed/biobanca/local/share/data/clinical_data_done.tsv"
msis <- read.table(msis_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
msis <- msis[,c(1, 23)]
colnames(msis) <- c("Model", "MSI_MSS")

res <- merge(merged, msis, by = "Model", all.x = TRUE)
res <- res[,c(1,2,3,4,11,7,8,9,5,6)]
rownames(res) <- res$Model

## sidedness 
res$Sidedness <- as.character(res$Sidedness)
res["CRC0277", "Sidedness"] <- "NONE"
for (i in seq(res$Model)) {
  if (res[i, "Sidedness"] == "SIGMOID COLON" | res[i, "Sidedness"] == "ANO" | 
      res[i, "Sidedness"] == "LEFT COLON" | res[i, "Sidedness"] == "RECTUM" |
      res[i, "Sidedness"] == "SPLENIC FLEXURE"){ 
    res[i, "Sidedness"] <- "LEFT COLON"
  } else if (res[i, "Sidedness"] == "CAECUM" | res[i, "Sidedness"] == "TRANSVERSE COLON" | 
             res[i, "Sidedness"] == "HEPATIC FLEXURE" | res[i, "Sidedness"] == "RIGHT COLON" ) {
    res[i, "Sidedness"] <- "RIGHT COLON"
  } else {
    res[i, "Sidedness"] <- NA
  }
}

## stage
res$Stage <- as.numeric(res$Stage)

## therapy before
## tutto ok caso che non ha info è il CRC0504

## MSI_MSS
is.na(res$MSI_MSS) <- NA
res$MSI_MSS <-  gsub('NT', NA, res$MSI_MSS)
res$MSI_MSS <- as.factor(res$MSI_MSS)

#KRAS #NRAS #BRAF
# tutto ok

#Sex
## già factor mancante CRC0277

# Age at collection
## mancanti CRC3034 e CRC0277
res$Age_at_Collection <- as.character(res$Age_at_Collection)
res$Age_at_Collection <- as.numeric(res$Age_at_Collection)

write.table(res, "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/chemiojul23_clinical_data_for_circos.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)
