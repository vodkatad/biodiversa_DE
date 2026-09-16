library(tidyverse)
library(ggplot2)
library(reshape2)
library(openxlsx)

tres <- read.table(file='/mnt/trcanmed/snaketree/stash/onf_etal_scores_all.tsv', sep="\t", quote="", header = TRUE, stringsAsFactors = FALSE)

tres$gen <- rownames(tres)
long <- melt(tres, id.vars='gen')
long$smodel <- substr(long$gen,0,7)
long$model <- substr(long$gen ,0,10)

#  no problem col .2 qui 
# tres[grepl('.2', tres$gen, fixed=T),]

long$type <- ifelse(grepl('LMX', long$model), 'Met', 'Pri')
ggplot(data=long[long$variable=='onf', ], aes(x=reorder(gen, value), y=value, fill=type))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                             axis.text.x=element_blank(),
                                                                                                             axis.ticks.x=element_blank())

long$type <- factor(long$type, levels=c('Pri', 'Met'))
ggplot(data=long, aes(x=type, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)

ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$type=='Pri', 'value'], myd[myd$type=='Met', 'value']))
}

## Media repliche e stat rifatto

meta_f <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meta_df<-read.table(meta_f, quote = "", sep = "\t", header = TRUE)
meta_df$sample_id_R <- gsub("-2", "", meta_df$sample_id_R)
meta_df$type <- gsub(".1", "", meta_df$type)
meta_df$type <- gsub(".2", "", meta_df$type)
 
meta <- meta_df %>% filter(type == "LMX_BASALE")

lmx_basale_ids <- meta$sample_id_R  

is_lmx <- grepl('LMX', tres$gen)
is_prx <- grepl('PRX', tres$gen)

# tieni: LMX che sono nella lista basale, OPPURE qualsiasi PRX
keep <- (is_lmx & tres$gen %in% lmx_basale_ids) | is_prx
cat("Totali:", nrow(tres), "\n")
cat("LMX_basali):", sum(is_lmx & tres$gen %in% lmx_basale_ids), "\n")
cat("LMX_non basali:", sum(is_lmx & !tres$gen %in% lmx_basale_ids), "\n")
cat("PRX:", sum(is_prx), "\n")
cat("Extraterrestri:", sum(!is_lmx & !is_prx), "\n")

tres <- tres[keep, ]

tres$smodel <- substr(tres$gen, 1, 10)
tres$type <- ifelse(grepl('LMX', tres$smodel), 'Met', 'Pri')

n_per_smodel <- table(tres$smodel)
reps <- n_per_smodel[n_per_smodel > 1]
print(reps)
print(table(tres[tres$smodel %in% names(reps), 'type']))

numeric_cols <- c('CBC','colon','esophagus','onf','stem','RSC','coreHRC','stem_index')
tres_mean <- tres %>%
  group_by(smodel, type) %>% ## qui type ridondante basta smodel perché CRCxxxxLMX o CRCxxxxPRX ma mi serve per dopo
  summarise(across(all_of(numeric_cols), mean), .groups = 'drop')

tres_mean <- as.data.frame(tres_mean)
write.table(tres_mean, "ssgsea_lmx_mediati.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
long <- melt(tres_mean, id.vars = c('smodel', 'type'), measure.vars = numeric_cols)
long$type <- factor(long$type, levels = c('Pri', 'Met'))

ggplot(data = long[long$variable == 'onf', ], aes(x = reorder(smodel, value), y = value, fill = type)) +
  geom_col() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

## riduzione a Onf stem_index coreHRC
longf <- long[long$variable %in% c('onf', 'coreHRC', 'stem_index'),]
# cor tests
ggplot(data = longf, aes(x = type, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  facet_wrap(~variable)



ggplot(data = long, aes(x = type, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.3) +
  facet_wrap(~variable)

# wilcoxon
ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s, ]
  cat(s, "\n")
  print(wilcox.test(myd[myd$type == 'Pri', 'value'], myd[myd$type == 'Met', 'value']))
}

# todo provare pairs (vedere bulk_moorman_sign_PriMets.R)
# TODO KRAS da Fra
# correlazioni tra coppie signature

ggplot(data=long, aes(x=reorder(smodel, value), y=value, fill=type))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.ticks.x=element_blank())+facet_wrap(~variable)

#valutazione KRAS
d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_update0625.tsv', sep='\t', header=T)
krasmut <- unique(d[d$KRAS!='wt', 'CASE'])

long$CASE <- substr(long$smodel, 1, 7)
long$KRAS <- ifelse(long$CASE %in% krasmut, 'Mut', 'WT')
long$KRAS <- factor(long$KRAS, levels=c('WT', 'Mut'))

ggplot(data=long[long$variable=='onf', ], aes(x=reorder(smodel, value), y=value, fill=KRAS))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                             axis.text.x=element_blank(),
                                                                                                             axis.ticks.x=element_blank())
# B1 violini
longf <- long[long$variable %in% c('onf', 'coreHRC', 'stem_index'),]
longf$variable <- as.character(longf$variable)
longf[longf$variable=='onf', 'variable'] <- 'OnF'
longf[longf$variable=='stem_index', 'variable'] <- 'stem index'
ggplot(data=longf, aes(x=KRAS, y=value))+ geom_violin(linewidth=0.5)+
  geom_boxplot(width=0.3, linewidth=0.5) + facet_wrap(~variable) + theme_minimal(base_size = 20)+
  ylab('ssGSEA')

ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$KRAS=='WT', 'value'], myd[myd$KRAS=='Mut', 'value']))
}


## valutazione per BRAF

d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_update0625.tsv', sep='\t', header=T)
brafmut <- unique(d[d$KRAS!='wt', 'CASE'])
brafmut <- brafmut[!is.na(brafmut)]

long$CASE <- substr(long$smodel, 1, 7)
long$BRAF <- ifelse(long$CASE %in% brafmut, 'Mut', 'WT')
long$BRAF <- factor(long$BRAF, levels=c('WT', 'Mut'))

table(long$BRAF)

long$BRAF <- factor(long$BRAF, levels=c('WT', 'Mut'))

ggplot(data=long[long$variable=='onf', ], aes(x=reorder(smodel, value), y=value, fill=BRAF))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                             axis.text.x=element_blank(),
                                                                                                             axis.ticks.x=element_blank())
ggplot(data=long, aes(x=BRAF, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)


ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$BRAF=='WT', 'value'], myd[myd$BRAF=='Mut', 'value']))
}

## PER DEG LMX

d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_update0625.tsv', sep='\t', header=T)
meta$CASE <- substr(meta$sample_id_R, 1, 7)
d$KRAS_status <- ifelse(d$KRAS == 'wt', 'WT', 'MUT')
res <- merge(meta, d[, c('CASE', 'KRAS_status')], by = 'CASE')
res <- res[,c("CASE", "sample_id_R", "KRAS_status")]
res <- res[!duplicated(res$sample_id_R),]
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
colnames(res) <- c("model", "KRAS")

sum(is.na(res$KRAS))
which(is.na(res$KRAS))
res[is.na(res$KRAS), ]
res <- res[!is.na(res$KRAS), ]

save.image(file = "plasticity_lmx.RData")


#write.table(res, "samples_data", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

## LMO

tres <- "/mnt/trcanmed/snaketree/stash/onf_etal_scores_all_pdo.tsv"
tres <- read.table(tres, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tres$gen <- rownames(tres)
tres$gen <- gsub(".2", "", tres$gen, fixed=T)

meta_f <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meta_df<-read.table(meta_f, quote = "", sep = "\t", header = TRUE)
meta_df$sample_id_R <- gsub("-2", "", meta_df$sample_id_R) # da fare anche in tres altrimenti sembra ci siano dei non basali
meta_df$type <- gsub(".1", "", meta_df$type, fixed=T)
meta_df$type <- gsub(".2", "", meta_df$type, fixed=T)

meta <- meta_df %>% filter(type == "LMO_BASALE")

lmo_basale_ids <- meta$sample_id_R  

is_lmo <- grepl('LMO', tres$gen)
is_prx <- grepl('PRX', tres$gen)

keep <- (is_lmo & tres$gen %in% lmo_basale_ids) | is_prx
cat("Totali:", nrow(tres), "\n")
cat("LMO_basali:", sum(is_lmo & tres$gen %in% lmo_basale_ids), "\n")
cat("LMO_non basali:", sum(is_lmo & !tres$gen %in% lmo_basale_ids), "\n")
cat("PRX:", sum(is_prx), "\n")
cat("Extraterrestri:", sum(!is_lmo & !is_prx), "\n")

tres <- tres[keep, ]

tres$smodel <- substr(tres$gen, 1, 10)
#tres$type <- ifelse(grepl('LMX', tres$smodel), 'Met', 'Pri')

n_per_smodel <- table(tres$smodel)
reps <- n_per_smodel[n_per_smodel > 1]
print(reps)
#print(table(tres[tres$smodel %in% names(reps), 'type']))

numeric_cols <- c('CBC','colon','esophagus','onf','stem','RSC','coreHRC','stem_index')
tres_mean <- tres %>%
  group_by(smodel) %>%
  summarise(across(all_of(numeric_cols), mean), .groups = 'drop')

tres_mean <- as.data.frame(tres_mean)
long <- melt(tres_mean, id.vars = c('smodel'), measure.vars = numeric_cols)
#long$type <- factor(long$type, levels = c('Pri', 'Met'))

#valutazione KRAS LMO
#d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/dtb_mutations_update0625.tsv', sep='\t', header=T)
krass <- "/mnt/trcanmed/snaketree/stash/pdo_xenturion_manual_tsg_info.xlsx"
krass <- read.xlsx(krass)
d <- krass[,c("CASE", "KRAS")]
krasmut <- unique(d[d$KRAS!='wt', 'CASE'])

long$CASE <- substr(long$smodel, 1, 7)
long$KRAS <- ifelse(long$CASE %in% krasmut, 'Mut', 'WT')
long$KRAS <- factor(long$KRAS, levels=c('WT', 'Mut'))

ggplot(data=long[long$variable=='onf', ], aes(x=reorder(smodel, value), y=value, fill=KRAS))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                              axis.text.x=element_blank(),
                                                                                                                axis.ticks.x=element_blank())

# B1 violini
longf <- long[long$variable %in% c('onf', 'coreHRC', 'stem_index'),]
longf$variable <- as.character(longf$variable)
longf[longf$variable=='onf', 'variable'] <- 'OnF'
longf[longf$variable=='stem_index', 'variable'] <- 'stem index'
ggplot(data=longf, aes(x=KRAS, y=value))+ geom_violin(linewidth=0.5)+
  geom_boxplot(width=0.3, linewidth=0.5) + facet_wrap(~variable) + theme_minimal(base_size = 20)+
  ylab('ssGSEA')

### waterfall PDOs per B1
lonf <- long[long$variable=='onf',]
sel <- c('CRC0099', 'CRC1634', 'CRC0188', 'CRC1979', 'CRC0116', 'CRC0095', 'CRC0123', 'CRC0169', 'CRC1241', 'CRC1566', 'CRC1477', 'CRC1239')

# unicode non si vede nel pannello plot ma nel pdf con cairo è ok
lonf$selected <- ifelse(lonf$CASE %in% sel, ifelse(lonf$value > 0, '↓','↑'), '')
lonf$vjust <- ifelse(lonf$value < 0, -0.05, 0.05)
g <- ggplot(data=lonf, aes(x=reorder(smodel, value), y=value, fill=KRAS))+geom_col()+
  scale_fill_manual(values=c('darkgrey', 'coral2'))+
  #scale_color_manual(values=c('darkgrey', 'black'))+
  geom_text(aes(label = selected, y = value + vjust), size=5)+
  theme_classic(base_size = 20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.line.x=element_blank())
 
ggsave(g, filename='~/whatever.pdf', device=cairo_pdf)

#####
ggplot(data=long, aes(x=KRAS, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)

ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$KRAS=='WT', 'value'], myd[myd$KRAS=='Mut', 'value']))
}


## valutazione per BRAF

brafs <- "/home/mferri/pdo_xenturion_manual_tsg_info.xlsx"
brafs <- read.xlsx(brafs)
d <- brafs[,c("CASE", "BRAF")]
brafmut <- unique(d[d$BRAF!='wt', 'CASE'])
brafmut <- brafmut[!is.na(brafmut)]

long$CASE <- substr(long$smodel, 1, 7)
long$BRAF <- ifelse(long$CASE %in% brafmut, 'Mut', 'WT')
long$BRAF <- factor(long$BRAF, levels=c('WT', 'Mut'))

table(long$BRAF)

long$BRAF <- factor(long$BRAF, levels=c('WT', 'Mut'))

ggplot(data=long[long$variable=='onf', ], aes(x=reorder(smodel, value), y=value, fill=BRAF))+geom_col() + theme(axis.title.x=element_blank(),
                                                                                                                axis.text.x=element_blank(),
                                                                                                                axis.ticks.x=element_blank())
ggplot(data=long, aes(x=BRAF, y=value))+ geom_violin()+
  geom_boxplot(width=0.3) + facet_wrap(~variable)

ss <- unique(long$variable)
for (s in ss) {
  myd <- long[long$variable == s,]
  print(s)
  print(wilcox.test(myd[myd$BRAF=='WT', 'value'], myd[myd$BRAF=='Mut', 'value']))
}


## PER DEG LMO 

krass <- "/home/mferri/pdo_xenturion_manual_tsg_info.xlsx"
krass <- read.xlsx(krass)
krass <- krass[,c("CASE", "KRAS")]
krass$KRAS <- ifelse(krass$KRAS == 'wt', 'WT', 'MUT')
meta_lmo <- meta_df %>% filter(type == "LMO_BASALE")
meta_lmo$CASE <- substr(meta_lmo$sample_id_R, 1, 7)
res <- merge(meta_lmo, krass, by = 'CASE')
res <- res[,c("CASE", "sample_id_R", "KRAS")]
res <- res[!duplicated(res$sample_id_R),]
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
colnames(res) <- c("model", "KRAS")
sum(is.na(res$KRAS))
which(is.na(res$KRAS))
res[is.na(res$KRAS), ]
res <- res[!is.na(res$KRAS), ]

## tolgo dai samples data questi che non sono nel sequenziamento di biodiversa per biobanca
del <- c("CRC0123LMO0A04015001R01000", "CRC0150LMO0A04014001R01000", "CRC0186LMO0B04024001R01000",
"CRC1078LMO0A05021001R01000", "CRC1078LMO0A05022001R01000", "CRC1090LMO0A05021001R01000",
"CRC1451LMO0A04008001R01000", "CRC1472LMO0A05013001R01000", "CRC1523LMO0A05021001R01000",
"CRC1523LMO0A05019001R01000", "CRC1620LMO0B04011001R01000")

res <- res %>% filter(!rownames(res) %in% del)

#write.table(res, "samples_data", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

