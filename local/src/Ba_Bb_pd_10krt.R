### selection of Ba_Bb_PD for keratin high and low

## carico meda e vsd per avere gli lmx e l'espressione
meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))
meda_f <- filter(meda_f, grepl("LMX_BASALE", type))

vsd <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'

## carico il tsv delle cheratine significative e prendo le 10 true
krt <- "/home/mferri/significant_krt.tsv"
krt <- read.table(krt, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
krtba <- (krt$genes)[krt$sign_Ac.vs.Ba == TRUE]
krtbb <- (krt$genes)[krt$sign_Ac.vs.Bb == TRUE]
krt_all <- unique(c(krtba, krtbb))

## carico connector e recist
connector <- "/scratch/trcanmed/connector/local/share/data/tsne_cetuxi_3w.tsv"
connector <- read.table(connector, quote = "", sep = "\t", header = T, stringsAsFactors = F)

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
names(recist)[names(recist) == "case"] <- "ShortID"

## merge tra connector e recist per filtrare Ac e PD
merged <- merge(connector, recist, by = "ShortID")
connector2 <- as.data.frame(cbind(merged$ShortID, merged$col, merged$classification))
colnames(connector2) <- c("case", "col", "classification")
con <- c("Ba", "Bb")
connector2 <- connector2 %>% filter(col %in% con)
connector2 <- connector2 %>% filter(classification == "PD")
connector_cases <- connector2$case

#filtro vsd per le cheratine gli lmx basali e le cheratine di interesse
vsd_subset <- vsd
vsd_subset <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_R, "symbol")]
vsd_subset <- vsd_subset[vsd_subset$symbol %in% krt_all,]
rownames(vsd_subset) <- vsd_subset$symbol
vsd_subset$symbol <- NULL
vsd_subset <- as.data.frame(t(vsd_subset))
vsd_subset$case <- substr(rownames(vsd_subset), 1, 7)
vsd_subset$genealogy <- rownames(vsd_subset)

##faccio il merge con il file di connector per ottenere solo i casi che mi interessano
##decido di mantenere il genealogy così può servire in futuro per dare a fra il campione giusto
merged2 <- merge(vsd_subset, connector2, by = "case")
setdiff(connector2$case, unique(merged2$case)) 
rownames(merged2) <- merged2$genealogy
merged2$genealogy <- NULL
merged2$case <- NULL
merged2$col <- NULL
merged2$classification <- NULL

write.table(merged2, file = "/home/mferri/Ba_Bb_PD_10_krt.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

## preparo un df per i quartili e calcolo i quartili
samples_krt <- merged2
quartiles_krt <- data.frame(matrix(ncol = 3, nrow = 10))
colnames(quartiles_krt) <- c("genes", "firstQ_25%", "thirdQ_75%")
quartiles_krt$genes <- colnames(samples_krt)

for (i in colnames(samples_krt)) {
  quartiles_krt[quartiles_krt$genes == i, "firstQ_25%"] <- quantile(samples_krt[,i], probs = c(0.25))
  quartiles_krt[quartiles_krt$genes == i, "thirdQ_75%"] <- quantile(samples_krt[,i], probs = c(0.75))
}
rownames(quartiles_krt) <- quartiles_krt$genes
#quartiles_krt$genes <- NULL
#quartiles_krt <- as.data.frame(t(quartiles_krt))
#prova <- as.data.frame(quantile(samples_krt$KRT80, probs = c(0.25,0.75)))

# ## preparo un df che sia identico a sample_krt solo vuoto per inserire se il valore exp
# ## è low del primo Q o high del terzo Q
# final <- data.frame(matrix(ncol = 10, nrow = 18))
# rownames(final) <- rownames(samples_krt)
# colnames(final) <- colnames(samples_krt)
# #x <- factor(c("low", "high"), levels = c("low", "high"))
# 
# for (i in rownames(samples_krt)) {
#   for (j in colnames(samples_krt)) {
#     for (k in rownames(quartiles_krt)) {
#       if(samples_krt[i, j] < quartiles_krt[k, "firstQ_25%"]) {
#         final[rownames(final) == i, colnames(final) == j] <- "low"
#       } else if (samples_krt[i, j] > quartiles_krt[k, "thirdQ_75%"]) {
#         final[rownames(final) == i, colnames(final) == j] <- "high"
#       }
#     }
#   }
# }
# EG con la struttura a tre cicli il problema è solo che manca il match delle keratine, nel senso che facendo un loop su tutti i quartili in questo momento
# noi confrontiamo i valori di expr della KRT80 al primo giro del for a riga 97 correttamente con la 80, ma dopo lo facciamo
# con i threshold della 7 e via dicendo...e alla fine nella matrice final metti sempre il confronto con i threshold dell'ultima krt del df quartiles_krt
# e non con quella giusta. Il fix è in realtà semplice perchè hai impostato bene sia i due for più esterni che le strutture dati, al posto del terzo
# for innestato faccio il controllo secco solo 'per la cheratina giusta'.
# Il sentore che ci fosse un problema l'ho avuto vedendo che erano tutte high o low: mi aspetto che ci siano molti valori di expr inclusi tra i due terzili
# (in effetti la maggioranza).
final_eg <- data.frame(matrix(ncol = 10, nrow = 9))
rownames(final_eg) <- rownames(samples_krt)
colnames(final_eg) <- colnames(samples_krt)
#x <- factor(c("low", "high"), levels = c("low", "high"))

for (i in rownames(samples_krt)) {
  for (j in colnames(samples_krt)) {
    wanted_krt <- j # variabile di appiggiosolo per chiarezza del codice
    first_q_25 <- quartiles_krt[rownames(quartiles_krt)==wanted_krt, 'firstQ_25%']
    third_q_75 <- quartiles_krt[rownames(quartiles_krt)==wanted_krt, 'thirdQ_75%']
    if (samples_krt[i, j] < first_q_25) {
      final_eg[rownames(final_eg) == i, colnames(final_eg) == j] <- "low"
    } else if (samples_krt[i, j] > third_q_75) {
      final_eg[rownames(final_eg) == i, colnames(final_eg) == j] <- "high"
    } else {
      final_eg[rownames(final_eg) == i, colnames(final_eg) == j] <- "middle"
    }
  }
}
final2 <- final_eg
#final2 <- as.data.frame(t(final2))
# if(samples_krt["CRC0080LMX0B02204TUMR05000", "KRT80"] < quartiles_krt["KRT80", "firstQ_25%"]) {
#   final["CRC0080LMX0B02204TUMR05000", "KRT80"] <- "low"
# } else if (samples_krt["CRC0080LMX0B02204TUMR05000", "KRT80"] > quartiles_krt["KRT80", "thirdQ_75%"]) {
#   final["CRC0080LMX0B02204TUMR05000", "KRT80"] <- "high"
# }
# prova <- data.frame(matrix(ncol = 10, nrow = 18))
# rownames(prova) <- rownames(samples_krt)
# colnames(prova) <- colnames(samples_krt)
# 
# if(samples_krt["CRC0420LMX0A02001TUMR06000", "KRT6A"] < quartiles_krt["KRT6A", "firstQ_25%"]) {
#   prova["CRC0420LMX0A02001TUMR06000", "KRT6A"] <- "low"
#   } else if (samples_krt["CRC0420LMX0A02001TUMR06000", "KRT6A"] > quartiles_krt["KRT6A", "thirdQ_75%"]) {
#    prova["CRC0420LMX0A02001TUMR06000", "KRT6A"] <- "high"
#   }
# 
# if(samples_krt["CRC0080LMX0B02204TUMR05000", "KRT6B"] < quartiles_krt["KRT6B", "firstQ_25%"]) {
#   prova["CRC0080LMX0B02204TUMR05000", "KRT6B"] <- "low"
#    } else if (samples_krt["CRC0080LMX0B02204TUMR05000", "KRT6B"] > quartiles_krt["KRT6B", "thirdQ_75%"]) {
#   prova["CRC0080LMX0B02204TUMR05000", "KRT6B"] <- "high"
#  }

# sum(!duplicated(final_binary[,"CRC0481LMX0B02001TUMR01000"]))==1
# sum(!duplicated(final_binary[,"CRC0080LMX0B02204TUMR05000"]))==1
# 
# length(unique(final_binary[,"CRC0481LMX0B02001TUMR01000"]))==1
# length(unique(final_binary[,"CRC0080LMX0B02204TUMR05000"]))==1

## faccio infiniti == e AND per vedere chi ha solo high o low

# for (k in rownames(final2)) {
#   if (final2[k, "KRT80"] == final2[k, "KRT7"] & final2[k, "KRT7"] == final2[k, "KRT86"] &
#     final2[k, "KRT86"] == final2[k, "KRT81"] & final2[k, "KRT86"] == final2[k, "KRT81"] & 
#     final2[k, "KRT86"] == final2[k, "KRT81"] & final2[k, "KRT81"] == final2[k, "KRT83"] &
#     final2[k, "KRT83"] == final2[k, "KRT6B"] & final2[k, "KRT6B"] == final2[k, "KRT6A"] &
#     final2[k, "KRT6A"] == final2[k, "KRT74"] & final2[k, "KRT74"] == final2[k, "KRT79"] &
#     final2[k, "KRT79"] == final2[k, "KRT16"]) {
#     final2[k, "check"] <- "chose"
#     } else {
#       final2[k, "check"] <- "lose"
#   }
# }

# il codice è corretto, io per non impazzire in caso di cambio di keratina ti propongo questo codice più generale con il fantomatico apply
check_all_high_low <- function(info) {
  if (all(info == "high")) {
    return("KRT_high")
  } else if (all(info == 'low')) {
    return("KRT_low")
  } else {
    return('KRT_ignavo')
  }
}

classification <- apply(final2, 1, check_all_high_low)

write.table(final2, file = "/home/mferri/Ac_Pd_high_low_krt.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

