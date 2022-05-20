library(ggplot2)
library(reshape)

pop_f <- snakemake@input[['pop']]
popc_f <- snakemake@input[['popc']]
lmx_f <- snakemake@input[['lmx']]
lmo_f <- snakemake@input[['lmo']]
lmh_f <- snakemake@input[['lmh']]
keep_f <- snakemake@input[['keepTRUE']]

plot_f <- snakemake@output[['plot']]
log_f <- snakemake@log[[1]]


load(snakemake@input[['Rimage']])

save.image('pippo.Rdata')
pop_freqs <- read.table(pop_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE, row.names=1)
pop_n <- read.table(popc_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE, row.names=1)

lmx_assign <- read.table(lmx_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)
lmo_assign <- read.table(lmo_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)
lmh_assign <- read.table(lmh_f, sep="\t", stringsAsFactors = FALSE,  header=TRUE)
keep_df <- read.table(keep_f, sep="\t", stringsAsFactors = FALSE, header= TRUE)
keep_model <- keep_df[keep_df$buoni == TRUE,'smodel'] 
#setwd('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')
#lmx="cris_vsd_lmx_nc_smodel.tsv" 
#lmo="cris_vsd_lmo_nc_smodel.tsv"
#lmh="cris_vsd_prediction_result_nc_lmh.tsv"
#pop_freqs <- read.table('freqs_cris.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
#pop_n <- read.table('counts_cris.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
#lmx_assign <- read.table(lmx, sep="\t", stringsAsFactors = FALSE, header=TRUE)
#lmo_assign <- read.table(lmo, sep="\t", stringsAsFactors = FALSE, header=TRUE)
#lmh_class <- read.table(lmh, sep="\t", stringsAsFactors = FALSE, header=TRUE)

#lmh_assign <- lmh_class[,c(1,2)]
#colnames(lmh_assign) <- c('genealogy', 'cris')

get_freq <- function(df, name, models, freq=TRUE) {
  df <- df[df$genealogy %in% models, , drop=FALSE]
  freqs <- as.data.frame(table(df$prediction))
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
  # since we are now filtering on subset of samples we could end up having no het even for X/O, so we change the check
  # that was adding  0 only to LMH
  #if (name == "LMH") { 
  # res$HET <- 0
  #}
  cols <- unique(c(colnames(res), 'NC', 'HET'))
  if (!'HET' %in% colnames(res)) {
    res$HET <- 0
  }
  if (!'NC' %in% colnames(res)) {
    res$NC <- 0
  }
  res <- res[,cols]
  return(res)
}

us <- rbind(get_freq(lmx_assign, 'LMX', keep_model), get_freq(lmo_assign, 'LMO', keep_model), get_freq(lmh_assign, 'LMH', keep_model))

usc <- rbind(get_freq(lmx_assign, 'LMX', keep_model, FALSE), get_freq(lmo_assign, 'LMO', keep_model, FALSE), get_freq(lmh_assign, 'LMH', keep_model, FALSE))


colnames(pop_freqs)[length(colnames(pop_freqs))] <- 'NC'
pop_freqs$HET <- rep(0, nrow(pop_freqs))
#colnames(pop_freqs) <- gsub('NS', 'NC', colnames(pop_freqs), fixed=TRUE)

pdata <- rbind(us, pop_freqs)
pdata$id <- rownames(pdata)
longpdata <- melt(pdata)
longpdata$variable <- factor(longpdata$variable, levels=c('CMS1', 'CMS2', 'CMS3', 'CMS4', 'NC', 'HET'))

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
ggsave(plot_f, h=4, w=4, units="in")



colnames(pop_n)[length(colnames(pop_n))] <- 'NC'
pop_n$HET <- rep(0, nrow(pop_n))
colnames(pop_n) <- gsub('CRIS.', 'CRIS-', colnames(pop_n), fixed=TRUE)
pdata <- rbind(usc, pop_n)
# We remove the HET class which we cannot have by definition in public datasets
pdata$HET <- NULL
# We decided to sum all the public dataset in a single one before going on with the chisq.
# We remove the first n rows from pdata, where n==nrow(usc), so the code will be more general and
# won't need many modifications to work on xeno only for magnum.

# Its a bit tortuous since we had them split, then put together, remove HET, and split again. 
# The code started with putting together to do the overall chisq, so I decided to keep it.
first_keep <- nrow(usc) + 1
if (first_keep <= nrow(pdata)) { # always better to be safe with seq, considering what seq(1,0) does
  public_nohet <- pdata[seq(first_keep, nrow(pdata)),]
  us_nohet <- pdata[seq(1, nrow(usc)),]
} else{
  save.image('pippo.Rdata')
  stop('You did not give me any info on public classification, something is wrong, check pippo')
}
sink(log_f)
print('Overall chisq')
chisq.test(pdata)
sink()


psum <- as.data.frame(apply(public_nohet, 2, sum))
#oursum <- as.data.frame(apply(us_nohet, 2, sum)) 
oursum <- t(us_nohet[rownames(us_nohet)=="LMX", , drop=FALSE]) # we keep only xeno as reference for chisq instead of summing everyone
chi_n <- cbind(psum, oursum)
sink(log_f, append=TRUE)
print('Summing pops chisq')
chisq.test(chi_n)
sink()


save.image('pippo.Rdata')