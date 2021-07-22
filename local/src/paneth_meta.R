
tmm_f <- snakemake@input[['tmm']]
output <- snakemake@output[['paneth_meta']]

tmm_f <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm_H.tsv.gz'
#wanted <- c('ATOH1','DEFA5','DEFA6','DLL1','GFI1','AREG', 'EREG','EGF','HBEGF','TGFA','BTC')

#wanted2 <- c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4","HER2","HER3","BTC")
#ATOH1, GFI1, SOX9,
#XBP1, DEFA5, DEFA6, LYZ, SPINK4, DLL1,
#and DLL4)

wanted <- c("ATOH1","GFI1","SOX9","XBP1","DEFA5","DEFA6","LYZ","SPINK4","DLL1","DLL4")

tmm <- read.table(gzfile(tmm_f), sep="\t", header=TRUE)

tmm <- tmm[wanted, ]

meta <- colMeans(log(tmm+1))

write.table(as.data.frame(meta), '/home/egrassi/panethmeta_biod5starok_selected.tsv', sep="\t",  quote=F)