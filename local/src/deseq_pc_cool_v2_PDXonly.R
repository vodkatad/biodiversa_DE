### need deseq_pc_cool_v2.R for all the functions etc

load("/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/pc_manual_plots.Rdata")
pca <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/type_PCx1_PCy2_PDX.pdf"
number1 <- as.numeric(1)
number2 <- as.numeric(2)
what <- "type"

vsd_PDX <- vsd[,grepl("X",colnames(vsd))]

p <- pc_cool(object = vsd_PDX, intgroup = what, pc1 = number1, pc2 = number2, plotfile=pca, legend=TRUE)

p

### leuco, using grid_pca_leuco.R
load("/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/leuco_score.Rdata")
df <- df[!grepl("H",df$sample),]
# df <- df[grepl("X",df$sample),]
df <- df[grepl("LMX_BASALE",df$type),]

color <- "Leucocyte"
### here launch the if part in grid_pca_leuco.R
p7
ggsave(filename="/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/type_PCx1_PCy2_LMXbasali-LeucoScore.pdf", plot=p7, dpi=300, height=120, width=180, units='mm')
getwd()

### check filter
filter <- read.table("/scratch/trcanmed/DE_RNASeq/local/share/data/magnum/DrugResponse_LMXfirslevel_trainTest.csv", header=TRUE)

table(filter$ircc_id_short %in% substr(df$sample,1,7))
# FALSE  TRUE 
# 3   228
notin <- filter[filter$ircc_id_short %in% substr(df$sample,1,7)==FALSE,"ircc_id_short"]
# CRC0680 CRC0018 CRC0246