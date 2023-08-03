library(data.table)
library(ggplot2)
library(pheatmap)
input<-snakemake@input[['data']]
input_s<-snakemake@input[['annot']]
data<-read.table(input,sep=',',header=TRUE,row.names = 1)
sign<-read.table(input_s,sep=',',header=TRUE,row.names=1)
print(nrow(data))
#row.names(sign)<-row.names(data)
#colnames(sign)<-colnames(data)
if(nrow(data)==0){
    prtin('empty')
  pdf(snakemake@output[['plot']])
  graphics.off()
}else{
pdf(snakemake@output[['plot']],width=12,height=12)
rg <- max(abs(data),na.rm=TRUE);
#pdf(snakemake@output[['plot']],width=12,height=12)
pheatmap(data,fontsize_row = 8,fontsize_col = 10,display_numbers = sign,nas_col='black',cluster_rows=FALSE,cluster_cols=FALSE,breaks = seq(-rg, rg, length.out = 100))
graphics.off()
}