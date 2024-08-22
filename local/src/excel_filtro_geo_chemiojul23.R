## excel per filtro geo

library(openxlsx)
library(tidyverse)

w3_f <- snakemake@input[["week3"]]
w6_f <- snakemake@input[["week6"]]
m29_f <- snakemake@input[["magnifici"]]
methy_f <- snakemake@input[["metilazione"]]
result <- snakemake@output[["tsv"]]

#w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
w3 <- read.table(w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

for (i in rownames(w3)) {
  if (w3[i, "type"] == "responder_1Q") {
    w3[i,"type"] <- 'sensitive tertile'
  } else {
    w3[i,"type"] <- 'resistant tertile'
  }
}

w3$genealogy <- rownames(w3)
w3 <- w3[,c(4,2)]

#w6_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23_6w/samples_data"
w6 <- read.table(w6_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

for (i in rownames(w6)) {
  if (w6[i, "type"] == "responder_1Q") {
    w6[i,"type"] <- 'sensitive tertile'
  } else {
    w6[i,"type"] <- 'resistant tertile'
  }
}

w6$genealogy <- rownames(w6)
w6 <- w6[,c(4,2)]

#m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/samples_data"
m29 <- read.table(m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

m29$genealogy <- rownames(m29)
m29 <- m29[,c(4,3)]

#methy_f <- "/mnt/trcanmed/snaketree/prj/methylomic/dataset/HR/list/list-selected_sample_info.tsv"
methy <- read.table(methy_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
methy <- methy[,c(1), drop(FALSE)]
colnames(methy) <- c("genealogy")

sgen <- w3
sgen$genealogy <- substr(sgen$genealogy, 1, 10)

methy <- merge(methy, sgen, by="genealogy", all.x=TRUE)
methy <- methy %>% filter(!duplicated(methy$genealogy))

# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "3 weeks of treatment")
addWorksheet(OUT, "6 weeks of treatment")
addWorksheet(OUT, "selected samples 6w of treat")
addWorksheet(OUT, "Methylation")

# Write the data to the sheets
writeData(OUT, sheet = "3 weeks of treatment", x = w3)
writeData(OUT, sheet = "6 weeks of treatment", x = w6)
writeData(OUT, sheet = "selected samples 6w of treat", x = m29)
writeData(OUT, sheet = "Methylation", x = methy)

# Export the file
saveWorkbook(OUT, result)
