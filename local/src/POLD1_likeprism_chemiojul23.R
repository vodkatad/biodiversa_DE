size <- 8
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- size* death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize, family='Arial'),
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


dds <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/dds.Rdata"
load(dds)
what <- "type"

data<-plotCounts(dds, "H_POLD1", intgroup=what, returnData=T)
data <- data[order(data$count),]
e2 <- ggplot(data, aes_string(x = what, y = "count"))
e3 <- e2 + geom_jitter(aes_string(shape = what, color = what),   position = position_jitter(0.2),size = 3) + stat_summary( aes_string(color = what), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.4, color="darkgreen")+theme_bw()+scale_y_continuous(trans='log10')+labs(color = "Xeno", x="Xeno", shape="Xeno", y="Log10(nreads)")+ggtitle("POLD1")+unmute_theme
ggsave(paste0("POLD1.eps"), width = 150, height = 100, units = "mm")
