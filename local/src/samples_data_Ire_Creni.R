library(stringr)
#ss<-read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/general/samples_data_general',header=TRUE,row.names = 1)
#ss<-ss[ss$model=='CRC0322',]
#prova<-ss[ss$trattamento == 'Creni_5uM_72h' | ss$trattamento == 'EGF0.1', ]
#meta<-read.table(meta_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE) 

meta_f <- snakemake@input[["meta"]]
ss<-read.table(meta_f,header=TRUE,row.names = 1,stringsAsFactors = FALSE)
res <- snakemake@output[["tsv"]]

path_ko<-'/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/ko_atoh1/general/samples_data'
ss_ko<-read.table(path_ko,header=TRUE,row.names = 1,stringsAsFactors = FALSE)
#ss_ko <- ss_ko[ !grepl("NO", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("0069", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("1139", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("1620", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("1502", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("542", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("CRC0322_CL7", row.names(ss_ko)),]
ss_ko <- ss_ko[!grepl("CRC0322_CL11", row.names(ss_ko)),]

colnames(ss_ko)<-colnames(ss)
ss<-rbind(ss,ss_ko)
#ss <-read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Ire_Creni/samples_data',header = TRUE,row.names = 1,stringsAsFactors = FALSE)

ss$trattamento<-ifelse(ss$trattamento == "EGF0.1", "EGF", ss$trattamento)
ss$trattamento<-ifelse(ss$trattamento == "CTX", "Cetux_72h", ss$trattamento)
ss$rep<-ifelse(ss$rep == "1", "WT_2", ss$rep)
ss$rep<-ifelse(ss$rep == "2", "WT_2", ss$rep)
ss$rep<-ifelse(ss$rep == "3", "WT_2", ss$rep)
ss$rep<-ifelse(ss$rep == "CL6", "KO_ATOH1", ss$rep)
ss$rep<-ifelse(ss$rep == "CL8", "KO_ATOH1", ss$rep)
ss$rep<-ifelse(ss$rep == "CL11", "KO_ATOH1", ss$rep)
ss$rep<-ifelse(ss$rep == "CL9", "KO_ATOH1", ss$rep)
ss$rep<-ifelse(ss$rep == "CAS9_2", "CAS9", ss$rep)
ss$rep<-ifelse(ss$rep == "CAS9_1", "CAS9", ss$rep)







#meta_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/ko_atoh1/ss_ko_atoh1.txt"
#ss <- read.table(meta_f,header=TRUE,row.names = 1)

#ss<-ss[ss$model=='CRC0327',]
#ss<-ss[ss$trattamento == 'Combo_72h' | ss$trattamento == 'Cetux_72h', ]
write.table(ss, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
