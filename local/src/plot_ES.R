library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(RColorBrewer)


rdata <- snakemake@input[["imagine"]]
data_f <- snakemake@input[["data"]]
out<-snakemake@output[['out']]

#rdata<-'/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/GSEA.Rdata'
load(rdata)
print(head(em))

p<-gseaplot2(em, geneSetID = "KEGG_NOTCH_SIGNALING_PATHWAY", title = "Enrichment Plot - KEGG_NOTCH_SIGNALING")
ggsave(p, file=out)
