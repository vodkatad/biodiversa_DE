#!/usr/bin/env Rscript

library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
genen <- args[[2]]
output1 <- args[[3]]
output2 <- args[[4]]

expr <- read.table(input, sep="\t", header=FALSE)
colnames(expr) <- c('lsample','gene')
expr$sample <- as.factor(substr(expr$lsample,0,7))
f <- as.data.frame(sapply(levels(expr$sample), function(x) { mean(expr[expr$sample==x,'gene']) }))
colnames(f) <- 'expr'
f$model <- rownames(f)

tt <- table(expr$sample)
dup <- names(tt[tt>1])
fdup <- expr[expr$sample %in% dup,]
fdup <- fdup[order(fdup[,'gene']),]
p <- ggplot(fdup, aes(x=reorder(sample, gene), y=gene)) +  geom_point() +theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(genen)
ggsave(output2, plot=p)

textSize <- 1
current_theme <-
  theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 12),
    #axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0,
                               size = 9),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y.left = unit(3,'mm'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(),
    legend.position = "top",
    legend.justification = "right",
    #legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(size = rel(textSize), face = "bold"),
    legend.text = element_text(size = rel(textSize)),
    legend.background = element_rect(size=0.5, linetype="solid", color="black"),
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    panel.border = element_blank()
  )
f$l <- log2(f$expr+1)
me <- median(f$l)
save.image('pippo.Rdata')
#p <- ggplot(f, aes(y=l,x=reorder(model, -l)))+geom_col(colour='black')+ylab('expr')+xlab("")+current_theme+scale_fill_manual(values=c("darkgrey","orange"))+ggtitle(genen)+geom_hline(aes(yintercept=me, linetype="1"), size=1,color="darkblue")+scale_x_discrete(expand=expansion(add=4))
p <- ggplot(f, aes(y=l,x=reorder(model, -l)))+geom_col()+ylab('expr')+xlab("")+current_theme+ggtitle(genen)+geom_hline(aes(yintercept=me, linetype="1"), size=1,color="darkblue")+scale_x_discrete(expand=expansion(add=4))
p <- p+scale_linetype_manual(name = "median", labels = "", values="solid") 

ggsave(output1, plot=p)

