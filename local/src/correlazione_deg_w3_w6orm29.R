library(ggplot2)
library(ggrastr)

w3_f <- snakemake@input[["week3"]]
w6_f <- snakemake@input[["week6"]]
m29_f <- snakemake@input[["magn29"]]
result <- snakemake@output[["jpg"]]
result_eps <- snakemake@output[["eps"]]
m29_result <- snakemake@output[["m29_jpg"]]
m29_result_eps <- snakemake@output[["m29_eps"]]
log_f <- snakemake@log[['log']]
  
#load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
#showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
#showtext_auto(enable = TRUE)

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

#w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
w3 <- read.table(w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#w6_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23_6w/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
w6 <- read.table(w6_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

w <- merge(w3, w6, by="row.names")

#nrow(m) # for the caption, need to put in a log file
ci <- cor.test(w$log2FoldChange.x, w$log2FoldChange.y) # idem
sink(log_f, append=TRUE)
"Correlation logFC w3 vs week6"
ci
sink()

w$padj_x <- ifelse(w$padj.x < 0.65, 0.65, w$padj.x)

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

y_breaks <- guess_ticks(w$log2FoldChange.x, fixed_min=-3, fixed_max=3)
x_breaks <- guess_ticks(w$log2FoldChange.y, fixed_min=-3, fixed_max=3)

p <- ggplot(data=w, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
  geom_point(size=0.1)+theme_bw()+xlab('logFC week 3')+ylab('logFC week 6')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+  
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
ggsave(p, file=result, width=4, height=4, units="in")

p <- ggplot(data=w, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
  geom_point(size=0.1) +theme_bw()+xlab('logFC week 3')+ylab('logFC week 6')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))
ggsave(p, file=result_eps, width=4, height=4, units="in")

## vs magnifici 29
#m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
m29 <- read.table(m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

w <- merge(w3, m29, by="row.names")

#nrow(m) # for the caption, need to put in a log file
ci <- cor.test(w$log2FoldChange.x, w$log2FoldChange.y) # idem
sink(log_f, append=TRUE)
"Correlation logFC w3 vs magnifici 29"
ci
sink()

w$padj_x <- ifelse(w$padj.x < 0.65, 0.65, w$padj.x)

y_breaks <- guess_ticks(w$log2FoldChange.x, fixed_min=-3, fixed_max=3)
x_breaks <- guess_ticks(w$log2FoldChange.y, fixed_min=-3, fixed_max=3)

p <- ggplot(data=w, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
  geom_point(size=0.1)+theme_bw()+xlab('logFC week 3')+ylab('logFC magnifici 29')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+  
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
ggsave(p, file=m29_result, width=4, height=4, units="in")

p <- ggplot(data=w, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj_x)))+
  geom_point(size=0.1) +theme_bw()+xlab('logFC week 3')+ylab('logFC magnifici 29')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))
ggsave(p, file=m29_result_eps, width=4, height=4, units="in")