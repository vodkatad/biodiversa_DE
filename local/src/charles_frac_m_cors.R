cr <- read.table('rnaseq_quality_control_information.tsv', header=T, quote="", comment.char="", sep="\t")

reads <- read.table(gzfile('/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/CharlesRiver/reads_info.tsv.gz'), sep="\t", header=TRUE)

cr$sample <- paste0(cr$Tumor_Designation, "_", cr$Model_Number, "_", cr$Sample_type)
m <- merge(reads, cr, by="sample")

stopifnot(nrow(m) == nrow(reads))

m$Alignment_rate <- as.numeric(gsub('%','', m$Alignment_rate))/100
m$Exonic_Reads <- as.numeric(gsub('%','', m$Exonic_Reads))/100

m$Mouse_Reads <- as.numeric(gsub('%','', m$Mouse_Reads))/100
m$us_m <- m$m_tot / m$Assigned
ggplot(data=m, aes(x=frac_assigned, y=Exonic_Reads, color=Strand))+geom_point()+theme_bw()+xlab('featureCounts assign bit')+ylab("% exonic cr")+theme(
  text = element_text(size = 20))+ggtitle("Pipeline QC - align/assign")

ggplot(data=m, aes(x=us_m, y=Mouse_Reads))+geom_point()+theme_bw()+xlab('% murine bit')+ylab("% murine cr")+theme(
  text = element_text(size = 20))+ggtitle("Pipeline QC - get ride of that mouse!")

d <- read.table('/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/metadata_charlesriver.tsv', sep="\t", header=TRUE)
mm <- merge(m, d, by.x="sample", by.y="id")
stopifnot(nrow(mm)==nrow(m))
library(ggsignif)
ggplot(data=mm, aes(x=Strand, y=cetuxi_perc_3w))+geom_boxplot(outlier.shape=NULL)+geom_jitter()
  +geom_signif(comparisons = list(c("firststrand", "unstranded")))+
  +theme(text = element_text(size = 20))