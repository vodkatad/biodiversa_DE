## file excel per i raw data dei volumi di cetux
cet <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3_buoni.tsv"
cet <- read.table(cet, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colnames(cet) <- c("Case", "Percentage volume cetuximab 3weeks")

write.xlsx(cet, file = "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3_buoni.xlsx")
