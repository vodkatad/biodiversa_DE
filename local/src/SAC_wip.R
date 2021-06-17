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

d<- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/merge_SAC-LMO_BASALE_tmm/matrix_tmm.tsv.gz'), sep="\t", header=T)
head(d)
dim(d)
rownames(d) <- d$X
d$X <- NULL
dd <- cor(d)
library(corrplot)
corrplot.mixed(dd)

decile <- function(x) {
  res <- cut(x, quantile(x, prob = seq(0, 1, length = 11), type = 5), include.lowest=TRUE )
  res <- as.factor(res)
  levels(res) <- seq(0,1, length=11)
  as.numeric(as.character(res))
}

deciles <- apply(d, 2, decile)
colnames(deciles) <- paste0(colnames(deciles),"_", "decile")

data <- cbind(d, deciles)
data$all_geq08 <- apply(deciles, 1, function(x) {all(x >= 0.8)} )
data$all_leq01 <- apply(deciles, 1, function(x) {all(x <= 0.1)} )
data <- data[order(-data$AURKA_decile),]
write.table(data,'SAC.tsv', sep="\t", quote=F)
data[data$all_geq08,]
data[data$all_leq01,]

data$avg_decile <- rowMeans(deciles)
data$model <- substr(rownames(data), 0, 7)

avg_avg_decile <- function(x, data) {
  myd <- data[data$model ==x,]
  return(mean(myd$avg_decile))
}

avgavg <- sapply(unique(data$model), avg_avg_decile, data)
davg <- as.data.frame(avgavg)
davg$model <- rownames(davg)
ctx <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w.tsv', sep="\t", header=F, row.names = NULL)
iri <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/may2020/irinotecan_w3.txt', sep="\t", header=T, row.names = NULL)
colnames(ctx) <- c('model', 'ctx')
colnames(iri) <- c('model', 'irino')

mctx <- merge(davg, ctx, by="model")
miri <- merge(davg, iri, by="model")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#https://slowkow.com/notes/ggplot2-color-by-density/
compare <- function(x, y,  log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  ggplot(d, aes(x=x, y=y)) +geom_point()+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+scale_color_manual(values=c('red', 'green'))+ggtitle(title)
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}

compare(mctx$avgavg, mctx$ctx, FALSE, 'avgDecile','cetuxi_vivo_3w','SAC')
compare(miri$avgavg, miri$irino, FALSE, 'avgDecile','irino_vivo_3w','SAC')
# only 90 vs cetuxi vivo is rather strange, no?
ctxpdo <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/cetuximab/pdo_cetuxi.tsv', sep="\t", header=T)
mctxo <- merge(davg, ctxpdo, by.x="model", by.y="case")
compare(mctxo$avgavg, mctxo$CTG_5000, FALSE, 'avgDecile','cetuxi_vitro_ATP5000','SAC')

###
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')
w <- c('CRC1629',
       'CRC1430',
       'CRC0080',
       'CRC0124',
       'CRC0186',
       'CRC1432')

load('LMO_BASALE-CEACAM5_tmm.png.Rdata')
wf <- function(data) {
f$track <- ifelse(f$model %in% w, 'yes','no')
f <- f[order(-f$l),]
f$lab <- f$model
f[f$track=="no",'lab'] <- ''
p <- ggplot(f, aes(y=l,x=reorder(model, -l), fill=track))+geom_col()+ylab('expr')+xlab("")+current_theme+ggtitle(genen)+geom_hline(aes(yintercept=me, linetype="1"), size=1,color="darkblue")+scale_fill_manual(values=c('gray','darkgoldenrod'))+scale_x_discrete(labels=f$lab)
p+scale_linetype_manual(name = "median", labels = "", values="solid") +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
load('LMO_BASALE-ERBB2_tmm.png.Rdata')
wf(f)
load('LMX_BASALE-ERBB2_tmm.png.Rdata')
wf(f)
load('LMX_BASALE-CEACAM5_tmm.png.Rdata')
wf(f)

########


d<- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/merge_SAC-LMX_BASALE_tmm/matrix_tmm.tsv.gz'), sep="\t", header=T)
head(d)
dim(d)
rownames(d) <- d$X
d$X <- NULL
dd <- cor(d)
library(corrplot)
corrplot.mixed(dd)

decile <- function(x) {
  res <- cut(x, quantile(x, prob = seq(0, 1, length = 11), type = 5), include.lowest=TRUE )
  res <- as.factor(res)
  levels(res) <- seq(0,1, length=11)
  as.numeric(as.character(res))
}

deciles <- apply(d, 2, decile)
colnames(deciles) <- paste0(colnames(deciles),"_", "decile")

data <- cbind(d, deciles)

w <- c('CRC0542',
'CRC0069',
'CRC1139',
'CRC0177',
'CRC1272',
'CRC0277',
'CRC1888',
'CRC0231'
)