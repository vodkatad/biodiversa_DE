data <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/MA_evil/vsd.tsv.gz', sep="\t", header=TRUE)
TOP <- 0.1

means <- apply(data, 1, mean)
med <- median(means)

### I obtain a data.table that contain only the genes  mean expression over the median
filtered <- as.data.frame(data[means > med,])

### I want to calculate the standard deviation
sds <- apply(filtered, 1, sd) 

### now we keep the top 10% variable genes 
sds <- sds[order(-sds)]
n <- length(sds)
keep <- head(sds, round(TOP*n))
keep_genes <- names(keep)
desd <- filtered[rownames(filtered) %in% keep_genes,]
ggplot(data=desd, aes(x=CRC1502_03_0, y=CRC1502_03_1_A))+geom_point()+geom_smooth(method='lm')

ggplot(data=desd, aes(x=CRC1502_09_0, y=CRC1502_09_1_C))+geom_point()+geom_smooth(method='lm')
cc <- cor(desd)
cor.test(desd$CRC1502_03_0, desd$CRC1502_03_1_A)
cor.test(desd$CRC1502_09_0, desd$CRC1502_09_1_C)

ggplot(data=desd, aes(x=CRC1502_09_0, y=CRC1502_09C_2_1))+geom_point()+geom_smooth(method='lm')
cor.test(desd$CRC1502_09_0, desd$CRC1502_09C_2_1)

jump_percentile <- function(d, s1, s2) {
  
  df <- data.frame(s1=d[,s1], s2=d[,s2], row.names=row.names(d))
  df$p1 <- cut(df$s1, quantile(df$s1, probs=seq(0, 1, by=0.1)), include.lowest=T, labels=F)
  df$p2 <- cut(df$s2, quantile(df$s2, probs=seq(0, 1, by=0.1)), include.lowest=T, labels=F)
  df$delta <- abs(df$p2-df$p1)
  #ggplot(data=df, aes(x=delta))+geom_histogram()+theme_bw(base_size=20)
  print(ggplot(data=df, aes(x=s1, y=s2, color=delta))+geom_point()+scale_color_viridis(discrete = FALSE)+theme_bw(base_size=20))
  return(nrow(df[df$delta>2,]))
}

jump_percentile(desd, 'CRC1502_03_0', 'CRC1502_03_1_A')
jump_percentile(desd, 'CRC1502_09_0', 'CRC1502_09C_2_1')

load('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/MA_evil/dds.Rdata')
library(DESeq2)
library(ggplot2)
df$time <- sapply(strsplit(rownames(df), "_"), function(x) {x[3]})
df$clone <- sapply(strsplit(rownames(df), "_"), function(x) {substr(x[2], 0, 2)})

pca_complicata <- function (matrix, ntop = 500){
  cv <- colVars(matrix)
  select <- order(cv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(cv)))]
  matrix <- matrix[, select]
  pca <- prcomp(matrix)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  return(list(pca = pca$x, pca_loading = pca$rotation, percentVar = percentVar))
}

tvsd <- t(assay(vsd))
list_pca <- pca_complicata(tvsd)

all(rownames(df)==rownames(list_pca$pca))
pd <- data.frame(PC1=list_pca$pca[,1],PC2=list_pca$pca[,2], time=df$time, clone=df$clone)

#pd$clone[is.na(pd$clone)] <- 'b'
#pd$time[is.na(pd$time)] <- '-1'

ggplot(data=pd, aes(x=PC1, y=PC2, color=time, shape=as.factor(clone)))+geom_point()+theme_bw(base_size=20)
