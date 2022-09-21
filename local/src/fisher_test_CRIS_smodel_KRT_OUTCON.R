## fisher test for CRIS* and krt high low

#connector <- "/scratch/trcanmed/connector/local/share/data/tsne_cetuxi_3w.tsv"
#connector <- read.table(connector, quote = "", sep = "\t", header = T, stringsAsFactors = F)

recist_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3.tsv"
recist <- read.table(recist_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
recist$classification <- NA

for (i in seq(length(recist$case))) {
  if (recist[i, "perc"] < -50.0){
    recist[i, "classification"] <- "OR"
  } else if (recist[i, "perc"] > -50.0 & recist[i, "perc"] < 35.0) {
    recist[i, "classification"] <- "SD"
  } else {
    recist[i, "classification"] <- "PD"
  }
}
names(recist)[names(recist) == "case"] <- "model"
# 
# merged <- merge(connector, recist, by  = "ShortID")
# merged <- as.data.frame(cbind(merged$ShortID, merged$col, merged$classification))
# colnames(merged) <- c("cases", "col", "classification")

krt <- "/home/mferri/OUTCON_Pd_high_low_krt.tsv"
krt <- read.table(krt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#krt$sample.names <- rownames(krt)
krt$model <- substr(rownames(krt), 1, 7)
krt$genealogy <- rownames(krt)

merged2 <- merge(recist, krt, by = "model")

cris <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/cris_vsd_ok_prediction_result_nc.tsv"
cris <- read.table(cris, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(cris)[names(cris) == "sample.names"] <- "genealogy"

merged <- merge(merged2, cris, by = "genealogy")
res <- as.data.frame(cbind(merged$model, merged$predict.label2, merged$KRT_classification))
colnames(res) <- c("model", "CRIS", "KRT_classification")
res$model <- as.character(res$model)
#res$col <- as.character(res$col)
res$CRIS <- as.character(res$CRIS)
res$KRT_classification <- as.character(res$KRT_classification)
res <- res[!duplicated(res), ]

write.table(res, file = "/scratch/trcanmed/DE_RNASeq/dataset/cris_deg_krt_connector/krt_cris_smodel_collapse.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# letter <- res
# crisa <- letter %>% filter(letter$CRIS == get("CRIS", envir = .env))
# other <- letter %>% filter(!letter$CRIS == get("CRIS", envir = .env))
# altibassi <- crisa
# high <- altibassi %>% filter(altibassi$KRT_classification == "KRT_high")
# low <- altibassi %>% filter(altibassi$KRT_classification == "KRT_low")
# altibassio <- other
# higho <- altibassio %>% filter(altibassio$KRT_classification == "KRT_high")
# lowo <- altibassio %>% filter(altibassio$KRT_classification == "KRT_low")
# ares <- data.frame(matrix(ncol = 2, nrow = 2))
# colnames(ares) <- c("High", "Low")
# rownames(ares) <- c("CRIS-A", "Other")
# ares[1,1] <- length(high$KRT_classification)
# ares[1,2] <- length(low$KRT_classification)
# ares[2,1] <- length(higho$KRT_classification)
# ares[2,2] <- length(lowo$KRT_classification)

get_tables_fisher <- function(CRIS) {
  letter <- res
  cris <- letter %>% filter(letter$CRIS == get("CRIS", envir = .env))
  other <- letter %>% filter(!letter$CRIS == get("CRIS", envir = .env))
  altibassi <- cris
  high <- altibassi %>% filter(altibassi$KRT_classification == "KRT_high")
  low <- altibassi %>% filter(!altibassi$KRT_classification == "KRT_high")
  altibassio <- other
  higho <- altibassio %>% filter(altibassio$KRT_classification == "KRT_high")
  lowo <- altibassio %>% filter(!altibassio$KRT_classification == "KRT_high")
  results <- data.frame(matrix(ncol = 2, nrow = 2))
  colnames(results) <- c("High", "Low")
  rownames(results) <- c("CRIS", "Other")
  results[1,1] <- length(high$KRT_classification)
  results[1,2] <- length(low$KRT_classification)
  results[2,1] <- length(higho$KRT_classification)
  results[2,2] <- length(lowo$KRT_classification)
  ress <- fisher.test(results)
  return(ress)
}

a <- get_tables_fisher("CRIS-A")
b <- get_tables_fisher("CRIS-B")
c <- get_tables_fisher("CRIS-C")
d <- get_tables_fisher("CRIS-D")
e <- get_tables_fisher("CRIS-E")
