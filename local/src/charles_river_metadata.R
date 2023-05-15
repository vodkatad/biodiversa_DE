library(reshape2)
setwd('~/work/def_targeted_sanger672_iorio/charlesriver/')


ss <- read.table('ss_charles.txt', sep="\t", header=TRUE)
soc <- read.table('soc.txt', sep="\t", header=TRUE, stringsAsFactors = FALSE)

rnaseq <- read.table('meta.txt', sep="\t", header=TRUE)

ss$id <- paste0(ss$Tumor.Designation, "_", ss$Tumor.Number, "_PDX")
m <- merge(rnaseq, ss[,-c(1,2)], by="id")
stopifnot(nrow(m) == nrow(rnaseq))

soc$id <- paste0(soc$Tissue.Type, "_", soc$Model.Number,  "_PDX")
soc <- soc[soc$id %in% m$id,]
soc <- soc[soc$Drug == "Cetuximab",]
stopifnot(length(soc$id), length(unique(soc$id)))

cetuxiv3w <- data.frame(id=soc$id, cetuxi_perc_3w=(soc$rel..Tumor.vol. -100))

# check vs approved dataset https://docs.google.com/presentation/d/17-_2z_vxQkoyKWIEOBUfEI5bQUyP4g5GqRWW4SijJG8/edit#slide=id.gc670f0b084_0_10
library(ggplot2)
ggplot(data=cetuxiv3w, aes(x=reorder(id, -perc_3w), y=perc_3w))+geom_bar(stat="identity")


mdata <- merge(m, cetuxiv3w, by='id', all.x=TRUE)
stopifnot(nrow(mdata) == nrow(rnaseq))
mdata$cetuxi_recist <- ifelse(mdata$cetuxi_perc_3w < -50, 'OR', ifelse(mdata$cetuxi_perc_3w > 35, 'PD','SD'))

ggplot(data=mdata[!is.na(mdata$cetuxi_perc_3w),], aes(x=reorder(id, -cetuxi_perc_3w), y=cetuxi_perc_3w, fill=cetuxi_recist))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=1))

ggplot(data=mdata[!is.na(mdata$cetuxi_perc_3w),], aes(x=reorder(id, -cetuxi_perc_3w), y=cetuxi_perc_3w, fill=cetuxi_recist))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=1))


# library(pheatmap)

# d <- read.table("~/work/def_targeted_sanger672_iorio/charlesriver/soc.txt", sep='\t', header=T, stringsAsFactors = F)
# d <- d[!is.na(d$WES),] 
# length(unique(d$Model.Number))
# length(unique(d$Drug))
# three <- d[d$WES & d$RNASEQ & d$SNP6,]
# two <- d[d$WES  & d$RNASEQ,]
# length(unique(three$Drug))
# length(unique(two$Drug))
# length(unique(two$Model.Number))
# length(unique(three$Model.Number))
# two$meas <- as.character(two$Measurement.day)
# #casted <- dcast(two, formula="Model.Number~Drug", fill="NA", value.var="meas", fun.aggregate=function(x){paste0(x, collapse=",")})
# casted <- dcast(two, formula="Model.Number~Drug", fill=NA_real_, value.var="rel..Tumor.vol.", fun.aggregate=function(x) {min(x)})
# rownames(casted) <- casted$Model.Number
# casted$Model.Number <- NULL
# dim(casted)
# cc <- two[,c('Model.Number','Tissue.Type')]
# cc <- unique(cc)
# cols <- data.frame(row.names=cc$Model.Number, Type=cc$Tissue.Type)
# pheatmap(t(casted), cluster_rows = FALSE, cluster_cols = F, annotation_col  = cols)
# ccc <- as.data.frame(table(two$Drug))
# rows <- data.frame(row.names=ccc$Var1, n=ccc$Freq)
# pheatmap(t(casted), cluster_rows = FALSE, cluster_cols = F, annotation_col  = cols, annotation_row=rows)
# hist(rows$n)
# 
# rownames(casted) <- casted$Model.Number
# casted$Model.Number <- NULL
# casted <- casted -1
# 
# y = v3/v0
# x = v3-v0/v0
# did I mean - 100?



write.table(mdata, file="metadata_charlesriver.tsv", sep="\t", quote=FALSE, row.names = FALSE)

## godot 15 05 23

#- model 233 has been treated for 314 days with Bevacizumab? Do you by any chance have also the value at 3 weeks? The number in this field is wrong, it should say 21 not 314!
  
#- model 243 have two different values for Irinotecan hydrochloride trihydrate, are they replicates? yes

#- the values are the ratio of the tumor volume at week 3 divided by week 0 (before treatment)? I know that you told Livio that they are relative to the starting point but  I wanted to be sure about the formula (since we use v3/v0 -1 usually and therefore have also negative values). Yes that is correct. We take day 0 as 100%

setwd('/mnt/trcanmed/snaketree/prj/magnum/local/share/data/charlesriver')

soc <- read.table('soc.txt', sep="\t", header=TRUE, stringsAsFactors = FALSE)
soc$id <- paste0(soc$Tissue.Type, "_", soc$Model.Number,  "_PDX")
soc <- soc[soc$Drug == "Irinotecan hydrochloride trihydrate",] # Oxaliplatin
soc$rel..Tumor.vol. <- as.numeric(soc$rel..Tumor.vol.)


nid <- as.data.frame(table(soc$id))
uniqueid <- nid[nid$Freq == 1, 'Var1']
avg243 <- mean(soc[soc$id=="CXF_243_PDX", 'rel..Tumor.vol.']) - 100

soc <- soc[soc$id %in% uniqueid,]
drug3w <- data.frame(id=c(soc$id, 'CXF_243_PDX'), perc_3w=c((soc$rel..Tumor.vol. -100), avg243))

drug3w$cetuxi_recist <- ifelse(drug3w$perc_3w < -50, 'OR', ifelse(drug3w$perc_3w > 35, 'PD','SD'))

ggplot(data=drug3w[!is.na(drug3w$perc_3w),], aes(x=reorder(id, -perc_3w), y=perc_3w, fill=cetuxi_recist))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=1))

