#### remove KRAS mut all together CTX72h/cronici xeno

library(ggplot2)
data <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/CSC-scores.tsv', header=T, sep="\t")
metadata <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', header=T, sep="\t", stringsAsFactors = FALSE)
data <- t(data)
m <- merge(data, metadata, by.x="row.names", by.y='id')
m$index <- m$RSC - m$Lgr5
m$type <- sapply(m$type, function(x) {y<-strsplit(x, '.', fixed=T)[[1]][1]; return(y[1])})
m$model <- substr(m$Row.names,0,10)

muts <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_vlookupAndrea_expandedPRLM.tsv', header=T, sep="\t")
muts <- unique(muts[,c('CASE','KRAS',"NRAS",'BRAF','PIK3CA')])

xeno <- m[grepl('LMX_cetux', m$type, fixed=TRUE),]
xeno$smodel <- substr(xeno$Row.names, 0, 7)
xeno <- xeno[xeno$smodel %in% muts$CASE,]
mutscores <- merge(xeno, muts, by.x="smodel", by.y="CASE")

ctx <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w.tsv', sep="\t", header=FALSE)
colnames(ctx) <- c('smodel', 'perc_cetuxi')

data2 <- merge(mutscores, ctx, by="smodel")

data2$time <- sapply( strsplit(data2$type, "_"), function(x) {x[3]})
print(ggplot(data=data2, aes(x=as.factor(time),y=index))+current_theme+theme(axis.text.x = element_text(angle=90))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2)))

#mmm <- merge(deltas, unique(data[,c('smodel','time')]), by="smodel")
#print(ggplot(data=mmm, aes(x=as.factor(time),y=delta))+current_theme+theme(axis.text.x = element_text(angle=90))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2)))


plastics2 <- function(data) {
  
  data$me <- paste0(data$type,"_", data$smodel)
  f <- as.data.frame(t(sapply(unique(data$me), function(x) { y<-data[data$me==x,]; c(mean(y[,'index']), unique(y[,'smodel'])) })))
  colnames(f) <- c('index','smodel')
  f$treat <- ifelse(grepl('_NT', rownames(f)),'NT','CTX')
  f$index <- as.numeric(as.character(f$index))
  #wilcox.test(formula=as.formula("index~treat"), data=f)
  
  du <- as.data.frame(table(f$smodel))
  du <- du[du$Freq ==2, "Var1"]
  ff <- f[f$smodel %in% du,]
  n <- ff[ff$treat=="NT",]
  t <- ff[ff$treat=="CTX",]
  all(t$smodel==n$smodel)
  print(ggplot(data=ff, aes(x=as.factor(treat),y=index))+current_theme+theme(axis.text.x = element_text(angle=90))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2)))
  print(wilcox.test(t$index, n$index, paired=TRUE))
  
  print(ggplot(ff) +
          geom_boxplot(aes(x = treat, y = index, group = treat))+
          geom_point(aes(x = treat, y = index)) +
          geom_line(aes(x = treat, y = index, group = smodel))) +current_theme
  
  
  deltas <- data.frame(delta = t$index - n$index, smodel=t$smodel, placebo_index = n$index, cetuxi_index = t$index)
  deltas$sign <- sign(deltas$placebo_index) * sign(deltas$cetuxi_index)
  
  deltas$plastic <- ifelse(deltas$sign == -1 & deltas$cetuxi_index > 0, 'plasticToRSC', ifelse(deltas$sign == -1 & deltas$cetuxi_index < 0, 'plasticToCSC', ifelse(deltas$placebo_index < 0, 'staticCBC','staticRSC')))
  
  mdf <- melt(deltas[,c('smodel','placebo_index', 'cetuxi_index')])
  meme <- merge(mdf, deltas, by="smodel")
  meme <- meme[match(mdf$smodel, meme$smodel),]
  mdf$plastic <- meme$plastic
  
  print(ggplot(mdf) +
          geom_boxplot(aes(x = variable, y = value, group = variable))+
          geom_point(aes(x = variable, y = value)) +
          geom_line(aes(x  = variable, y = value, group = smodel)) +current_theme+facet_wrap(~plastic))
  
  m2 <- merge(mdf, unique(data[,c('smodel','cetuxi')]), by="smodel")
  m2$cetuxi <- m2$cetuxi * 100
  print(ggplot(data=m2, aes(x=as.factor(plastic),y=cetuxi))+current_theme+theme(axis.text.x = element_text(angle=90))+geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2)))
  print(wilcox.test(m2[m2$plastic == "plasticToRSC",'cetuxi'],m2[m2$plastic != "plasticToRSC",'cetuxi']))
  print(wilcox.test(m2[m2$plastic == "plasticToCSC",'cetuxi'],m2[m2$plastic != "plasticToCSC",'cetuxi']))
  # delta in response q[1] q[0]
  return(m2)
}



data2 <- data2[data2$KRAS=="wt" & data2$NRAS=="wt" & data2$BRAF=="wt" & data2$PIK3CA=="wt",]
dim(data2)
length(unique(data2$smodel))

m2 <- plastics2(data2)


###
m2$class <- ifelse(m2$cetuxi < -50, 'OR', ifelse(m2$cetuxi > 35, "PD", "SD")) #Recist 3w: OR -50% /  SD / +35% PD
table(m2$plastic, m2$class)
pd <- as.data.frame(table(m2$plastic, m2$class))
ggplot(pd, aes(x="", y=Freq, fill=Var2))+ geom_bar(width = 1, stat = "identity")+facet_wrap(~Var1)+current_theme+coord_polar("y", start=0)

#https://github.com/tidyverse/ggplot2/issues/2815
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE

ggplot(pd, aes(x="", y=Freq, fill=Var2))+ geom_bar(width = 1, stat = "identity")+facet_wrap(~Var1, scales = "free")+current_theme+cp


#ggplot(pd, aes(x="", y=Freq, fill=Var2))+ geom_bar(width = 1, stat = "identity")+facet_wrap(~Var1)+current_theme+  theme(aspect.ratio = 1)+coord_polar("y", start=0)#

m <- as.matrix(table(m2$plastic, m2$class))
chisq.test(m)
mm <- m[c(1,2),]
chisq.test(mm)

##########

# compute delta for each model
delta <- function(data) {
  
  data$me <- paste0(data$type,"_", data$smodel)
  f <- as.data.frame(t(sapply(unique(data$me), function(x) { y<-data[data$me==x,]; c(mean(y[,'index']), unique(y[,'smodel'])) })))
  colnames(f) <- c('index','smodel')
  f$treat <- ifelse(grepl('_NT', rownames(f)),'NT','CTX')
  f$index <- as.numeric(as.character(f$index))
  #wilcox.test(formula=as.formula("index~treat"), data=f)
  
  du <- as.data.frame(table(f$smodel))
  du <- du[du$Freq ==2, "Var1"]
  ff <- f[f$smodel %in% du,]
  n <- ff[ff$treat=="NT",]
  t <- ff[ff$treat=="CTX",]
  all(t$smodel==n$smodel)
  
  deltas <- data.frame(delta = t$index - n$index, smodel=t$smodel, placebo_index = n$index, cetuxi_index = t$index)
  return(deltas)
}

delta_index <- delta(data2)

# compute PIS for each model

d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)
#wanted <- data.frame(gs=c('ATOH1','DEFA5','DEFA6','DLL1','GFI1','AREG', 'EREG','EGF','HBEGF','TGFA','BTC'))
paneth <- data.frame(gs=c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4"))

colnames(paneth) <- 'gs'
paneth$hgs <- paste0('H_', paneth$gs)

fpkm <- d[rownames(d) %in% paneth$hgs,]
fpkm_pdo_treat <- fpkm[,colnames(fpkm) %in% data2$Row.names]
pseudoc <- 1
lfpkm_pdo_treat <- log(fpkm_pdo_treat+pseudoc)

ave <- colMeans(lfpkm_pdo_treat) # run this way, changes a lot
#ave <- colMeans(fpkm_pdo_treat) 
paneth_treat <- data.frame(ave =ave, id = names(ave)) 
paneth_treat_m <- merge(metadata, paneth_treat, by="id")
paneth_treat_m$model <- substr(paneth_treat_m$id, 0, 7)
paneth_treat_m <- paneth_treat_m[order(paneth_treat_m$model, paneth_treat_m$type),]

fc <- function(model, data) {
  d <- data[data$model == model,]
  if (nrow(d) == 4) {
    fc <- mean(c(d[1, 'ave'] / d[2, 'ave'], d[4, 'ave'] / d[3, 'ave'])) # due to order TODO
  } else {
    fc <- d[1, 'ave'] / d[2, 'ave']
  }
  return(fc)
}

panethIndScore <- sapply(unique(paneth_treat_m$model), fc, paneth_treat_m)

pis <- data.frame(smodel=names(panethIndScore), pis=panethIndScore)

mpis <- merge(pis, delta_index, by='smodel')

# TODO indagare
#> dim(mpis)
#[1] 51  5
#> length(unique(pis$smodel))
#[1] 52
#> length(unique(delta_index$smodel))
#[1] 51

# correlate

# how do we put here 'extreme responders'?
plot(mpis$pis, mpis$delta)

cor.test(mpis$pis, mpis$delta)

mpisc <- merge(mpis, ctx, by='smodel')

mpisc$class <- ifelse(mpisc$perc_cetuxi < -50, 'OR', ifelse(mpisc$perc_cetuxi > 35, 'PD' ,'SD') )

ggplot(data=mpisc, aes(x=pis,y=delta, color=class))+geom_point()+current_theme

mpisc$lcet <- log(mpisc$perc_cetuxi+100)
ggplot(data=mpisc, aes(x=pis,y=delta, color=lcet))+geom_point()+current_theme+scale_color_gradient(low = "yellow", high = "red", na.value = NA)


