library(tidyverse)
library(ggplot2)
library(ggrastr)

grigi <- snakemake@input[["grey"]]
m29 <- snakemake@input[["magnifici"]] 
plot <- snakemake@output[["pdf"]]
log_f <- snakemake@log[["log"]]

#load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
# showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
# showtext_auto(enable = TRUE)

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- size* death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize),#, family='Arial'),
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
#axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),


# function that given values to be plotted on an axis will return:
# vector of breaks, trying to guess which max will be the best one
# this will be used as scale_y_continuous(breaks=  and as ylim(min, max) to have the - also limits-c()
# last tick at the extremity of the axis.
# other parameter is n. of ticks
guess_ticks <- function(values, nticks=5, fixed_max=NULL, fixed_min=0) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- ceiling(vmax)
  } else {
    round_max <- fixed_max
  }
  my_breaks <- seq(fixed_min, round_max, length.out=nticks)
  return(my_breaks)
}


#m29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
m29 <- read.table(m29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
m29$genes <- gsub("H_", "", rownames(m29))
names(m29)[names(m29)=="log2FoldChange"] <- "log2FoldChange_m29"
m29 <- m29[,c("genes", "log2FoldChange_m29")]

#grigi <- "/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
grigi <- read.table(grigi, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
grigi$genes <- gsub("H_", "", rownames(grigi))
names(grigi)[names(grigi)=="log2FoldChange"] <- "log2FoldChange_grigi"
grigi <- grigi[,c("genes", "log2FoldChange_grigi")]

merged <- merge(grigi, m29, by="genes")

sink(log_f, append=TRUE)
"Correlation magnifici29 grey"
cor.test(merged$log2FoldChange_m29, merged$log2FoldChange_grigi, method = "spearman")
sink()

y_breaks <- guess_ticks(merged$log2FoldChange_grigi, fixed_min=-5, fixed_max=5)
x_breaks <- guess_ticks(merged$log2FoldChange_m29, fixed_min=-3, fixed_max=3)

p <- ggplot(merged, aes(x=log2FoldChange_m29, y=log2FoldChange_grigi)) + 
  rasterise(geom_point(size = 1), dpi=300)+
  geom_smooth(method=lm, se=FALSE)+ # ratio is ... 1 becomes 0.939
  unmute_theme+theme(legend.position="none") + xlab('log2FoldChange_m29')+ylab('log2FoldChange_grigi')+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

#ggsave(p, file="/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/scatter_m29_grigi.pdf", width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")
ggsave(p, file=plot, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")


