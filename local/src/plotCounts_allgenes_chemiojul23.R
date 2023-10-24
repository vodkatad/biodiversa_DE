library(tidyverse)
library(DESeq2)
library(ggplot2)
library(reshape)
library(showtext)

dds_f <- snakemake@input[["dds"]]
pvalues_f <- snakemake@input[["deseq"]]
results <- snakemake@output[["plotcounts"]]
results_mute <- snakemake@output[["plotcounts_mute"]]
results_tsv <- snakemake@output[["tsv"]]

#load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
showtext_auto(enable = TRUE)

#textSize <- size * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- s * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- (size) * death_conversion_dpi96

unmute_theme <- theme(
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size=0.508/0.564), # origin of this ratio is honestly not known, empirical
  axis.ticks = element_line(color = "black", size=0.508/0.564),
  axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),
  panel.background = element_blank()
)
# #axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),


# # function that given values to be plotted on an axis will return:
# # vector of breaks, trying to guess which max will be the best one
# # this will be used as scale_y_continuous(breaks=  and as ylim(min, max) to have the - also limits-c()
# # last tick at the extremity of the axis.
# # other parameter is n. of ticks
# guess_ticks <- function(values, nticks=5, fixed_max=NULL, fixed_min=0) {
#   vmax <- max(values)
#   if (is.null(fixed_max)) { 
#     round_max <- ceiling(vmax)
#   } else {
#     round_max <- fixed_max
#   }
#   my_breaks <- seq(fixed_min, round_max, length.out=nticks)
#   return(my_breaks)
# }

#dds_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/dds.Rdata"
dds_file <- dds_f
#load("/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/dds.Rdata")
load(dds_file)
geni <- c("FBXO18", "RING1", "SUMO1", "UBE2I", "BRCA1", "RRP1",
          "FBXO5", "RFWD3",
          "UCHL3", "PARPBP", "BLM")
geni <- paste0("H_", geni)


#d <- plotCounts(dds, geni, intgroup="type", returnData=TRUE)

df_list <- list()
for (i in geni) {
  df_list[[i]] <- plotCounts(dds, gene=i, intgroup="type", returnData=TRUE)
  colnames(df_list[[i]]) <- c(i, "type")
}

combined_df <- bind_cols(df_list)
combined_df$type...4 <- NULL
combined_df$type...6 <- NULL
combined_df$type...8 <- NULL
combined_df$type...10 <- NULL
combined_df$type...12 <- NULL
combined_df$type...14 <- NULL
combined_df$type...16 <- NULL
combined_df$type...18 <- NULL
combined_df$type...20 <- NULL
combined_df$type...22 <- NULL
names(combined_df)[names(combined_df)=="type...2"] <- "type"

reshaped <- melt(combined_df, id = c("type"))
reshaped$variable <- gsub("H_", "", reshaped$variable)

save.image("plot_grid.Rdata")

#y_breaks <- guess_ticks(data$`RING1_log2(TPM+1)`, fixed_min=-100, fixed_max=5000)
#x_breaks <- guess_ticks(data$`FBXO18_log2(TPM+1)`, fixed_min=1, fixed_max=8)
red <- rgb(165, 0, 25, max = 255)
blue <- rgb(30, 85, 130, max = 255)

r_s <- ggplot(reshaped, aes(x = variable, y = value, color = type)) + geom_point(position =position_dodge(width = 0.7))+
  unmute_theme + scale_color_manual(values=c(red, blue))+
  aes(x = fct_inorder(variable))+theme(axis.title.x = element_blank())
ggsave(r_s, filename = results, width=200*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

#print("change column name to log10(nreads)")

r_m <- ggplot(reshaped, aes(x = variable, y = value, color = type)) + geom_point(position =position_dodge(width = 0.7))+
  unmute_theme + scale_color_manual(values=c("#a50019", "#1e5582"))+
  aes(x = fct_inorder(variable))+ theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                                          axis.text.x=element_blank(), axis.text.y=element_blank())
ggsave(r_m, filename = results_mute, width=150*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

#pvalues_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
pvalues <- read.table(pvalues_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pvalues$genes <- rownames(pvalues)
pvalues <- pvalues %>% filter(genes %in% geni)
pvalues <- pvalues[match(geni, pvalues$genes), ]
pvalues$genes <- NULL
rownames(pvalues) <- gsub("H_", "", rownames(pvalues))
pvalues <- pvalues[,c(5,6)]

write.table(pvalues, file=results_tsv, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)