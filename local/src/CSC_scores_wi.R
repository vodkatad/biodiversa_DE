library(ggplot2)

data <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/CSC-scores.tsv', header=T, sep="\t")
metadata <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', header=T, sep="\t", stringsAsFactors = FALSE)

data <- t(data)
m <- merge(data, metadata, by.x="row.names", by.y='id')
m$index <- m$RSC - m$Lgr5

m$type <- sapply(m$type, function(x) {y<-strsplit(x, '.', fixed=T)[[1]][1]; return(y[1])})
m$model <- substr(m$Row.names,0,10)

#################33
th <- function() {
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
current_theme
}
current_theme <- th()
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
compare <- function(x, y, index, log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  med <- median(index)
  d$med <- ifelse(index > med, 'larger','smaller')
  ggplot(d, aes(x=x, y=y, color=med)) +geom_point()+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+scale_color_manual(values=c('red', 'green'))+ggtile(title)
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
mutscores$fwt <- mutscores$KRAS=='wt' & mutscores$BRAF=='wt' & mutscores$NRAS=='wt' & mutscores$PIK3CA=="wt"
library(ggsignif)
ggplot(data=mutscores, aes(x=KRAS_bin,y=index, fill=KRAS_bin))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("mut", "wt")))

ggplot(data=mutscores, aes(x=BRAF_bin,y=index, fill=BRAF_bin))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("mut", "wt")))


ggplot(data=mutscores, aes(x=fwt,y=index, fill=fwt))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index")+geom_signif(comparisons = list(c("TRUE", "FALSE")))+xlab("Quadruple wt")
###
me <- read.table('/mnt/trcanmed/snaketree/prj/pdx_methylation/dataset/v2/heatmap/cluval/k5_samples-clusters_division_switched.tsv', sep="\t", header=TRUE)

basali <- m[grepl('X_BASALE',m$type),]
basali$model <- substr(basali$Row.names,0,10)
mergeme <- merge(basali, me, by.x="model", by.y="row.names")
ggplot(data=mergeme, aes(x=as.factor(cluster),y=index, fill=as.factor(cluster)))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index vs Methylation clusters")+geom_signif(comparisons = list(c("4", "5")))+theme( legend.position = "None")


### cris
cris <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/cris_tmm_0.2_classes.tsv', header=T, sep="\t")

mecri <- merge(m, cris, by.x="Row.names", by.y="genealogy")
ggplot(data=mecri, aes(x=cris,y=index, fill=cris))+geom_boxplot()+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Index vs CRIS")+theme( legend.position = "None")

########## 

compare <- function(x, y, index, log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  med <- median(index)
  d$med <- ifelse(index > med, 'larger','smaller')
  print(ggplot(d, aes(x=x, y=y, color=med)) +geom_point(size=2)+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+scale_color_manual(values=c('red', 'green'))+ggtitle(title))
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}

call_compare <- function(x, data) {
  d <- data[data$type == x,]
  #x, y, log, nx, ny, title
  compare(x=d$Lgr5, y=d$RSC, index=d$index, log=FALSE, nx='CBC', ny='RSC', title=x)
}

call_compare('LMX_BASALE', m)
call_compare('LMH', m)
###
# LMX vs LMH

data <- m[m$type %in% c('LMH','LMX_BASALE'),]
data$smodel <- substr(data$model,0,7)
data$class <- substr(data$model,8,10)

xe <- data[data$class=="LMX",]
h <- data[data$class=="LMH",c('smodel','index')]


tt <- table(xe$smodel)
dup <- names(tt[tt>1])
fdup <- xe[xe$smodel %in% dup,]
fdup <- fdup[order(fdup[,'index']),]
p <- ggplot(fdup, aes(x=reorder(smodel, index), y=index)) +  geom_point() +theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Index vs xeno replicates")

f <- as.data.frame(sapply(unique(xe$smodel), function(x) { mean(xe[xe$smodel==x,'index']) }))
colnames(f) <- 'index'
f$smodel <- rownames(f)

xeh <- merge(f, h, by="smodel")
colnames(xeh)[2] <- "xeno"
colnames(xeh)[3] <- "human"
mxeh <- melt(xeh)
wilcox.test(xeh$xeno, xeh$human, paired=TRUE)
ggplot(data=mxeh, aes(x=variable,y=value, fill=variable))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2))+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Met X vs H")+theme( legend.position = "None")+ylab('index')+xlab("sample")

compare2 <- function(x, y, log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y)
  ggplot(d, aes(x=x, y=y)) +geom_point(size=2)+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle(title)+geom_smooth(method = lm)
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}

compare2(x=xeh$xeno, y=xeh$human, log=FALSE, nx="xeno", ny="human",title="Mets")
### PRX vs PRH

data <- m[m$type %in% c('PRH','PRX_BASALE'),]
data$smodel <- substr(data$model,0,7)
data$class <- substr(data$model,8,10)

xe <- data[data$class=="PRX",]
h <- data[data$class=="PRH",c('smodel','index')]


tt <- table(xe$smodel)
dup <- names(tt[tt>1])
fdup <- xe[xe$smodel %in% dup,]
fdup <- fdup[order(fdup[,'index']),]
p <- ggplot(fdup, aes(x=reorder(smodel, index), y=index)) +  geom_point() +theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle("Index vs xeno replicates")
print(p)
f <- as.data.frame(sapply(unique(xe$smodel), function(x) { mean(xe[xe$smodel==x,'index']) }))
colnames(f) <- 'index'
f$smodel <- rownames(f)

xeh <- merge(f, h, by="smodel")
colnames(xeh)[2] <- "xeno"
colnames(xeh)[3] <- "human"
mxeh <- melt(xeh)
wilcox.test(xeh$xeno, xeh$human, paired=TRUE)
ggplot(data=mxeh, aes(x=variable,y=value, fill=variable))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2))+current_theme+theme(axis.text.x = element_text(angle=90))+ggtitle("Pri X vs H")+theme( legend.position = "None")+ylab('index')+xlab("sample")

compare2(x=xeh$xeno, y=xeh$human, log=FALSE, nx="xeno", ny="human",title="Primary")


#### cetuxi

cetuxi <- m[grepl('LMX_cetux_72h', m$type),c('index','model','type')]
ggplot(data=cetuxi, aes(x=type,y=index, fill=type))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2))+current_theme+theme(axis.text.x = element_text(angle=90))+theme( legend.position = "None")+ylab('index')+xlab("Treatment")

cetuxi$smodel <- substr(cetuxi$model, 0,7)
cetuxi$me <- paste0(cetuxi$type,"_", cetuxi$smodel)

f <- as.data.frame(t(sapply(unique(cetuxi$me), function(x) { y<-cetuxi[cetuxi$me==x,]; c(mean(y[,'index']), unique(y[,'smodel'])) })))

colnames(f) <- c('index','smodel')
f$treat <- ifelse(grepl('_NT', rownames(f)),'NT','CTX')
f$index <- as.numeric(as.character(f$index))
wilcox.test(formula=as.formula("index~treat"), data=f)

du <- as.data.frame(table(f$smodel))
du <- du[du$Freq ==2, "Var1"]

ff <- f[f$smodel %in% du,]
n <- ff[ff$treat=="NT",]
t <- ff[ff$treat=="CTX",]
all(t$smodel==n$smodel)
wilcox.test(t$index, n$index, paired=TRUE)


