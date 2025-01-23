## analysis of genes present in m29_grigi and not in the original m29

m29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
m29 <- read.table(m29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tot_m29 <- m29
m29 <- m29 %>% filter(padj < 0.05)
m29 <- m29 %>% filter(abs(log2FoldChange) > 0.5849625)
m29$genes <- gsub("H_", "", rownames(m29))
names(m29)[names(m29)=="log2FoldChange"] <- "log2FoldChange_m29"
m29 <- m29[,c("genes", "log2FoldChange_m29")]

grigi <- "/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
grigi <- read.table(grigi, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tot_grigi <- grigi
grigi <- grigi %>% filter(padj < 0.05)
grigi <- grigi %>% filter(abs(log2FoldChange) > 0.5849625)
grigi$genes <- gsub("H_", "", rownames(grigi))
names(grigi)[names(grigi)=="log2FoldChange"] <- "log2FoldChange_grigi"
names(grigi)[names(grigi)=="padj"] <- "padj_grigi"
grigi <- grigi[,c("genes", "log2FoldChange_grigi", "padj_grigi")]
grigi <- grigi[order(grigi$padj_grigi),]

res <- merge(grigi, m29, by="genes", all.x=TRUE, all.y = TRUE)
res$log2FoldChange_grigi[is.na(res$log2FoldChange_grigi)] <- "not_differential_in_DEG_grigi"
res$log2FoldChange_m29[is.na(res$log2FoldChange_m29)] <- "not_differential_in_DEG_m29"
rownames(res) <- res$genes

only_m29 <- res %>% filter(!log2FoldChange_m29 == "not_differential_in_DEG_m29")
only_grigi <- res %>% filter(!log2FoldChange_grigi == "not_differential_in_DEG_grigi")

both <- res %>% filter(!log2FoldChange_grigi == "not_differential_in_DEG_grigi")
both <- both %>% filter(!log2FoldChange_m29 == "not_differential_in_DEG_m29")
both$log2FoldChange_grigi <- as.numeric(both$log2FoldChange_grigi)
both$log2FoldChange_m29 <- as.numeric(both$log2FoldChange_m29)

for (i in rownames(both)) {
  if (both[i,"log2FoldChange_grigi"]*both[i,"log2FoldChange_m29"]>0) {
    both[i,"concordanza"] <- "concordi"
  } else {
    both[i, "concordanza"] <- "discordi"
  }
}

tot_grigi$genes <- gsub("H_", "", rownames(tot_grigi))
names(tot_grigi)[names(tot_grigi)=="log2FoldChange"] <- "log2FoldChange_grigi_tot"
names(tot_grigi)[names(tot_grigi)=="padj"] <- "padj_grigi_tot" 
tot_grigi <- tot_grigi[,c("genes", "log2FoldChange_grigi_tot", "padj_grigi_tot")]
tot_m29$genes <- gsub("H_", "", rownames(tot_m29))
names(tot_m29)[names(tot_m29)=="log2FoldChange"] <- "log2FoldChange_m29_tot"
names(tot_m29)[names(tot_m29)=="padj"] <- "padj_m29_tot" 
tot_m29 <- tot_m29[,c("genes", "log2FoldChange_m29_tot", "padj_m29_tot")]

sign_grigi_wm29 <- merge(only_grigi, tot_m29, by="genes")

setdiff(only_grigi$genes, sign_grigi_wm29$genes)
#pvalue
# [1] "AC004837.2" "ATG9B"      "CAMK2B"     "CERS1"      "CYP26A1"    "FBLN1"      "FGF19"      "FN1"        "KRT5"       "MAP2"      
# [11] "MAPK8IP2"   "MORN3"      "NLRP1"      "NPEPL1"     "NPTXR"      "NTSR1"      "PAX5"       "PPBP"       "ROBO2"      "SLC16A6"   
# [21] "SMARCA1"    "TLL2"       "TM4SF4"     "TM6SF2"     "VGF"        "XPNPEP2" 
#padj
# "ATG9B" "CERS1" "VGF"

names(sign_grigi_wm29)[names(sign_grigi_wm29)=="log2FoldChange_m29"] <- "in_m29"
#sign_grigi_wm29$in_m29 <- "no_DEG"
rownames(sign_grigi_wm29) <- sign_grigi_wm29$genes

for (i in rownames(sign_grigi_wm29)) {
  if (sign_grigi_wm29[i, "padj_m29_tot"] < 0.05) {
    sign_grigi_wm29[i,"padj"] <- "sign"
  } else {
    sign_grigi_wm29[i,"padj"] <- "non_sign"
  }
}

for (i in rownames(sign_grigi_wm29)) {
  if (abs(sign_grigi_wm29[i, "log2FoldChange_m29_tot"]) > 0.5849625) {
    sign_grigi_wm29[i,"differenziale"] <- "diff"
  } else {
    sign_grigi_wm29[i,"differenziale"] <- "no_diff"
  }
}

sign_grigi_wm29$why <- paste0(sign_grigi_wm29$differenziale, "_", sign_grigi_wm29$padj)

# for (i in rownames(sign_grigi_wm29)) {
#   if (sign_grigi_wm29[i,"differenziale"]=="no" & sign_grigi_wm29[i,"padj"]=="non_significativo") {
#     sign_grigi_wm29[i,"why"] <- "no_diff_no_sign"
#   } else if (sign_grigi_wm29[i,"differenziale"]=="sì" & sign_grigi_wm29[i,"padj"]=="non_significativo") {
#     sign_grigi_wm29[i,"why"] <- "diff_no_sign"
#   } else if (sign_grigi_wm29[i,"differenziale"]=="sì" & sign_grigi_wm29[i,"padj"]=="significativo") {
#     sign_grigi_wm29[i,"why"] <- "diff_sign"
#   } else {
#     sign_grigi_wm29[i,"why"] <- "no_diff_sign"
#   }
# }

sign_grigi_wm29$padj <- NULL
sign_grigi_wm29$differenziale <- NULL

for (i in rownames(sign_grigi_wm29)) {
  if (sign_grigi_wm29[i, "in_m29"]=="not_differential_in_DEG_m29") {
    sign_grigi_wm29[i, "in_m29"] <- "no_DEG"
  } else {
    sign_grigi_wm29[i, "in_m29"] <- "DEG"
  }
}

sign_grigi_wm29 <- sign_grigi_wm29[order(sign_grigi_wm29$padj_grigi),]

write.xlsx(sign_grigi_wm29, file="/home/mferri/differenziali_grigi_in_m29.xlsx")

sign_m29_wgrigi <- merge(only_m29, tot_grigi, by="genes")

setdiff(only_m29$genes, sign_m29_wgrigi$genes)
# [1] "AC004837.2" "ATG9B"      "CAMK2B"     "CERS1"      "CYP26A1"    "FBLN1"      "FGF19"      "FN1"        "KRT5"       "MAP2"      
# [11] "MAPK8IP2"   "MORN3"      "NLRP1"      "NPEPL1"     "NPTXR"      "NTSR1"      "PAX5"       "PPBP"       "ROBO2"      "SLC16A6"   
# [21] "SMARCA1"    "TLL2"       "TM4SF4"     "TM6SF2"     "VGF"        "XPNPEP2" 

names(sign_m29_wgrigi)[names(sign_m29_wgrigi)=="log2FoldChange_grigi"] <- "in_grigi"
sign_m29_wgrigi$in_grigi <- "no_DEG"
rownames(sign_m29_wgrigi) <- sign_m29_wgrigi$genes

for (i in rownames(sign_m29_wgrigi)) {
  if (sign_m29_wgrigi[i, "padj_grigi_tot"] > 0.05) {
    sign_m29_wgrigi[i,"padj"] <- "non_significativo"
  } else {
    sign_m29_wgrigi[i,"padj"] <- "significativo"
  }
}

for (i in rownames(sign_m29_wgrigi)) {
  if (abs(sign_m29_wgrigi[i, "log2FoldChange_grigi_tot"] > 0.5849625)) {
    sign_m29_wgrigi[i,"differenziale"] <- "sì"
  } else {
    sign_m29_wgrigi[i,"differenziale"] <- "no"
  }
}

for (i in rownames(sign_m29_wgrigi)) {
  if (sign_m29_wgrigi[i,"differenziale"]=="no" & sign_m29_wgrigi[i,"padj"]=="non_significativo") {
    sign_m29_wgrigi[i,"why"] <- "no_diff_no_sign"
  } else if (sign_m29_wgrigi[i,"differenziale"]=="sì" & sign_m29_wgrigi[i,"padj"]=="non_significativo") {
    sign_m29_wgrigi[i,"why"] <- "diff_no_sign"
  } else {
    sign_m29_wgrigi[i,"why"] <- "no_diff_sign"
  }
}

sign_m29_wgrigi$padj <- NULL
sign_m29_wgrigi$differenziale <- NULL
