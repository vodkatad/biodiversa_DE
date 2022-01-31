library(ggplot2)
library(reshape)

pop_f <- snakemake@input[['pop']]
popc_f <- snakemake@input[['popc']]
lmx_f <- snakemake@input[['lmx']]
lmo_f <- snakemake@input[['lmo']]
lmh_f <- snakemake@input[['lmh']]
plot_f <- snakemake@output[['plot']]
log_f <- snakemake@log[[1]]


load(snakemake@input[['Rimage']])

save.image('pippo.Rdata')
pop_freqs <- read.table(pop_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE, row.names=1)
pop_n <- read.table(popc_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE, row.names=1)

lmx_assign <- read.table(lmx_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)
lmo_assign <- read.table(lmo_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)
lmh_class <- read.table(lmh_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)

#setwd('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')
#lmx="cris_vsd_lmx_nc_smodel.tsv" 
#lmo="cris_vsd_lmo_nc_smodel.tsv"
#lmh="cris_vsd_prediction_result_nc_lmh.tsv"
#pop_freqs <- read.table('freqs_cris.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
#pop_n <- read.table('counts_cris.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
#lmx_assign <- read.table(lmx, sep="\t", stringsAsFactors = FALSE, header=TRUE)
#lmo_assign <- read.table(lmo, sep="\t", stringsAsFactors = FALSE, header=TRUE)
#lmh_class <- read.table(lmh, sep="\t", stringsAsFactors = FALSE, header=TRUE)

lmh_assign <- lmh_class[,c(1,2)]
colnames(lmh_assign) <- c('genealogy', 'cris')

get_freq <- function(df, name, freq=TRUE) {
  freqs <- as.data.frame(table(df$cris))
  if (freq) {
  freqs$frac <- freqs$Freq / sum(freqs$Freq)
  } else {
    freqs$frac <- freqs$Freq
  }
  freqs$Freq <- NULL
  rownames(freqs) <- freqs$Var1
  freqs$Var1 <- NULL
  res <- t(freqs)
  rownames(res) <- name
  res <- as.data.frame(res)
  if (name == "LMH") {
    res$HET <- 0
  }
  return(res)
}

us <- rbind(get_freq(lmx_assign, 'LMX'), get_freq(lmo_assign, 'LMO'), get_freq(lmh_assign, 'LMH'))

usc <- rbind(get_freq(lmx_assign, 'LMX', FALSE), get_freq(lmo_assign, 'LMO', FALSE), get_freq(lmh_assign, 'LMH', FALSE))


colnames(pop_freqs)[length(colnames(pop_freqs))] <- 'NC'
pop_freqs$HET <- rep(0, nrow(pop_freqs))
colnames(pop_freqs) <- gsub('CRIS.', 'CRIS-', colnames(pop_freqs), fixed=TRUE)

pdata <- rbind(us, pop_freqs)
pdata$id <- rownames(pdata)
longpdata <- melt(pdata)
longpdata$variable <- factor(longpdata$variable, levels=c('CRIS-A', 'CRIS-B', 'CRIS-C', 'CRIS-D', 'CRIS-E', 'NC', 'HET'))

longpdata$order <- 0
longpdata[grepl('GSE', longpdata$id), 'order'] <- 1
longpdata[longpdata$id == "LMH", 'order'] <- 2
longpdata[longpdata$id == "LMX", 'order'] <- 3
longpdata[longpdata$id == "LMO", 'order'] <- 4
ggplot(data=longpdata, aes(x=reorder(id, order), y=value, fill=variable))+geom_col()+
unmute_theme +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
scale_fill_brewer(palette="Dark2")+
xlab('Dataset')+
ylab('Fraction of samples')
ggsave(plot_f)



colnames(pop_n)[length(colnames(pop_n))] <- 'NC'
pop_n$HET <- rep(0, nrow(pop_n))
colnames(pop_n) <- gsub('CRIS.', 'CRIS-', colnames(pop_n), fixed=TRUE)
pdata <- rbind(usc, pop_n)
pdata$HET <- NULL
sink(log_f)
chisq.test(pdata)
sink()