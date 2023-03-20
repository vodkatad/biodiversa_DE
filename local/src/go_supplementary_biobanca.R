## supplementary table go biobanca

xoup <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE_up.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
xodown <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMX_BASALE_down.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ohup <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMH_up.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ohdown <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMO_BASALE.vs.LMH_down.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
xhup <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMX_BASALE.vs.LMH_up.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
xhdown <- read.table("/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_type_cutoff0.05-LMX_BASALE.vs.LMH_down.tsv",
                   quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

write.xlsx(xoup, file="GO_results_supplementary.xlsx", sheetName="LMO.vs.LMX_up", row.names=FALSE)
write.xlsx(xodown, file="GO_results_supplementary.xlsx", sheetName="LMO.vs.LMX_down", append=TRUE, row.names=FALSE)
write.xlsx(ohup, file="GO_results_supplementary.xlsx", sheetName="LMO.vs.LMH_up", append=TRUE, row.names=FALSE)
write.xlsx(ohdown, file="GO_results_supplementary.xlsx", sheetName="LMO.vs.LMH_down", append=TRUE, row.names=FALSE)
write.xlsx(xhup, file="GO_results_supplementary.xlsx", sheetName="LMX.vs.LMH_up", append=TRUE, row.names=FALSE)
write.xlsx(xhdown, file="GO_results_supplementary.xlsx", sheetName="LMX.vs.LMH_down", append=TRUE, row.names=FALSE)


# Create a blank workbook
OUT <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(OUT, "PDXT.vs.PDX_up")
addWorksheet(OUT, "PDXT.vs.PDX_down")
addWorksheet(OUT, "PDXT.vs.LMH_up")
addWorksheet(OUT, "PDXT.vs.LMH_down")
addWorksheet(OUT, "PDX.vs.LMH_up")
addWorksheet(OUT, "PDX.vs.LMH_down")

# Write the data to the sheets
writeData(OUT, sheet = "PDXT.vs.PDX_up", x = xoup)
writeData(OUT, sheet = "PDXT.vs.PDX_down", x = xodown)
writeData(OUT, sheet = "PDXT.vs.LMH_up", x = ohup)
writeData(OUT, sheet = "PDXT.vs.LMH_down", x = ohdown)
writeData(OUT, sheet = "PDX.vs.LMH_up", x = xhup)
writeData(OUT, sheet = "PDX.vs.LMH_down", x = xhdown)

# Export the file
saveWorkbook(OUT, "/scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca/Clustering_Cutoff0.05_LFC0.584_ok/GO_results_supplementary.xlsx")
