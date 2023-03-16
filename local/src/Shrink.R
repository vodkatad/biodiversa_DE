input <- snakemake@input[['chemo']]
outputTsv <- snakemake@output[['f']]
outputPlot <- snakemake@output[['p']]

#read the file with the volumes of folfiri and cetuximab treated models after 3 weeks and 6 weeks 
fol_cet <- read.table(input, sep="\t", header=TRUE)
fol_cet

#extract the first 3 columns of fol_cet and create a new data frame for the folfiri treated models only
fol <- fol_cet[, c(1, 2, 3)]

#split the observations into 2 groups: fast shrinkage and slow shrinkage + exclude all the other observations

#delete all the observations with NA (missing values of volume at week 6)
fol <- na.omit(fol)


for(i in 1:nrow(fol)){
  if(fol[i, "folfiri.vol_3w"] < 0 && fol[i, "folfiri.vol_6w"] < 0){
    fol[i, "shrink"] <- "fast" 
  }
  else if(fol[i, "folfiri.vol_3w"] > 0 && fol[i, "folfiri.vol_6w"] < 35 && fol[i, "folfiri.vol_6w"] < fol[i, "folfiri.vol_3w"]){
    fol[i, "shrink"] <- "slow"
  }
  else{
    fol[i, "shrink"] <- "non-shrink"
  }
}

#plot a graph for the fast shrinkage, slow shrinkage, and non-shrinkage observations
# we need to add the starting 0 for clarity
library(reshape)
fol0 <- data.frame(model=fol$model, folfiri.vol_0w=rep(0, nrow(fol)))
fol_with0 <- merge(fol, fol0, by="model")                    
longfol <- melt(data=fol_with0, id.vars=c("model", "shrink"), measure.vars=c('folfiri.vol_3w','folfiri.vol_6w', 'folfiri.vol_0w'))

library(ggplot2)
# reorder the factor to have 0w before
longfol$variable <- factor(longfol$variable, levels=c('folfiri.vol_0w','folfiri.vol_3w', 'folfiri.vol_6w'))
ggplot(data=longfol, aes(y=value, x=variable, group=model))+geom_line(alpha=0.3, aes(col=shrink))+scale_color_manual(values=c('blue', 'green', 'red'))+geom_point()+theme_bw()
ggsave(outputPlot)

write.table(fol_with0, file=outputTsv, quote=FALSE, sep='\t')