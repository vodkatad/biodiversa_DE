library(ggplot2)

data <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/CSC-scores.tsv', header=T, sep="\t")
metadata <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', header=T, sep="\t", stringsAsFactors = FALSE)

data <- t(data)
m <- merge(data, metadata, by.x="row.names", by.y='id')
m$index <- m$RSC - m$Lgr5

m$type <- sapply(m$type, function(x) {y<-strsplit(x, '.', fixed=T)[[1]][1]; return(y[1])})

#################33
textSize <- 1.5
current_theme <-
  theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = rel(1.8)),
    axis.text.x = element_text(size=rel(1.7)),
    axis.text.y = element_text(angle = 0,
                               size = rel(1.7)),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y.left = unit(3,'mm'),
    
    legend.position = "top",
    legend.justification = "right",
    #legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(size = rel(textSize), face = "bold"),
    legend.text = element_text(size = rel(1.2)),
    legend.background = element_rect(size=0.5, linetype="solid", color="black"),
    plot.title = element_text(
      face = "bold",
      size = rel(2),
      hjust = 0.5
    ),
    panel.border = element_blank(),
    plot.caption = element_text(size=rel(1))
  )

###
ggplot(data=m, aes(x=RSC))+geom_histogram()+facet_wrap(~type)+current_theme

ggplot(data=m, aes(x=type,y=RSC))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("RSC")
ggplot(data=m, aes(x=type,y=Lgr5))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Lgr5")
ggplot(data=m, aes(x=type,y=index))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#https://slowkow.com/notes/ggplot2-color-by-density/
compare <- function(x, y, log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  ggplot(d, aes(x=x, y=y, color=density)) +geom_point()+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+scale_color_viridis_c()+ggtitle(title)
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}

pe <- cor.test(m$RSC, m$Lgr5)

ggplot(m, aes(x=RSC, y=Lgr5, color=type)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")

mm <- m[m$type=="LMX_BASALE" & !is.na(m$cetuxi),]

pe <- cor.test(mm$RSC, mm$cetuxi); ggplot(mm, aes(x=RSC, y=cetuxi)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")
pe <- cor.test(mm$Lgr5, mm$cetuxi); ggplot(mm, aes(x=Lgr5, y=cetuxi)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")
pe <- cor.test(mm$index, mm$cetuxi); ggplot(mm, aes(x=index, y=cetuxi)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")



mm <- m[m$type=="LMX_BASALE" & !is.na(m$irino),]

pe <- cor.test(mm$RSC, mm$irino); ggplot(mm, aes(x=RSC, y=irino)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")
pe <- cor.test(mm$Lgr5, mm$irino); ggplot(mm, aes(x=Lgr5, y=irino)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")
pe <- cor.test(mm$index, mm$irino); ggplot(mm, aes(x=index, y=irino)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")


#
muts <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_vlookupAndrea_expandedPRLM.tsv', header=T, sep="\t")
muts <- unique(muts[,c('CASE','KRAS',"NRAS",'BRAF','PIK3CA')])
xeno <- m[m$type == "LMX_BASALE",]
xeno$model <- substr(xeno$Row.names, 0, 7)
xeno <- xeno[xeno$model %in% muts$CASE,]
mutscores <- merge(xeno, muts, by.x="model", by.y="CASE")

o <- as.data.frame(table(mutscores$Row.names))
oo <- as.character(o[o$Freq>1,"Var1"])
mutscores[mutscores$Row.names %in% oo,]

mutscores <- mutscores[-c(190,192,446,448),] # probably best to keep M only...
o <- as.data.frame(table(mutscores$Row.names))
oo <- as.character(o[o$Freq>1,"Var1"])
mutscores[mutscores$Row.names %in% oo,]

pe <- cor.test(mutscores$index, mutscores$cetuxi); ggplot(mutscores, aes(x=index, y=cetuxi, color=KRAS)) +geom_point()+current_theme+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle("SC Scores")

mutscores$KRAS_bin <- ifelse(mutscores$KRAS == "wt", "wt", "mut")
mutscores$BRAF_bin <- ifelse(mutscores$BRAF == "wt", "wt", "mut")
mutscores$wt <- mutscores$KRAS=='wt' & mutscores$BRAF=='wt' & mutscores$NRAS=='wt' & mutscores$PIK3CA=="wt"
library(ggsignif)
ggplot(data=mutscores, aes(x=KRAS_bin,y=index))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("mut", "wt")))

ggplot(data=mutscores, aes(x=BRAF_bin,y=index))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("mut", "wt")))


ggplot(data=mutscores, aes(x=wt,y=index))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("TRUE", "FALSE")))
