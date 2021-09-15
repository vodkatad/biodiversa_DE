data <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/tmm.tsv.gz'), sep="\t", header=T)

metadata <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', header=T, sep="\t", stringsAsFactors = FALSE)



muts <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_vlookupAndrea_expandedPRLM.tsv', header=T, sep="\t")
muts <- unique(muts[,c('CASE','KRAS',"NRAS",'BRAF','PIK3CA')])

wanted <- metadata[grepl('LMX_BASALE', metadata$type, fixed=TRUE),]

data <- data[, colnames(data) %in% wanted$id]
# ctx <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w.tsv', sep="\t", header=FALSE)
# colnames(ctx) <- c('smodel', 'perc_cetuxi')

plot(log2(as.numeric(data[rownames(data)=='H_IGF2',])+1),log2(as.numeric(data[rownames(data)=='H_MIR483',])+1))
cor.test(log2(as.numeric(data[rownames(data)=='H_IGF2',])+1),log2(as.numeric(data[rownames(data)=='H_MIR483',])+1))


muts58 <- read.table("~/Book2.txt", sep="\t", header=T)
muts58$wt <- ifelse(muts58$muts == "Quadruple WT", '4wt','4mut')


pdata <- data.frame(IGF=log2(as.numeric(data[rownames(data)=='H_IGF2',])+1), MIR483=log2(as.numeric(data[rownames(data)=='H_MIR483',])+1), lmodel=colnames(data))

pdata$smodel <- substr(pdata$lmodel, 0,7)

mm <- merge(pdata, muts58, by="smodel")
#w <- read.table('/tmp/o')
#plot(log2(as.numeric(data2[rownames(data2)=='H_IGF2',])+1),log2(as.numeric(data2[rownames(data2)=='H_MIR483',])+1))
#cor.test(log2(as.numeric(data2[rownames(data2)=='H_IGF2',])+1),log2(as.numeric(data2[rownames(data2)=='H_MIR483',])+1))
#cor.test(log2(as.numeric(data[rownames(data)=='H_IGF2',])+1),log2(as.numeric(data[rownames(data)=='H_MIR483',])+1))
plot(log2(as.numeric(data[rownames(data)=='H_IGF2',])+1),log2(as.numeric(data[rownames(data)=='H_MIR483',])+1))

ggplot(data=mm, aes(x=IGF,y=MIR483, color=wt))+geom_point()+current_theme

cor.test(mm$IGF, mm$MIR483)