
library(dplyr)

ctx <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w.tsv', sep="\t", header=F, row.names = NULL)
iri <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/may2020/irinotecan_w3.txt', sep="\t", header=T, row.names = NULL)
colnames(ctx) <- c('model', 'ctx')
colnames(iri) <- c('model', 'irino')
km <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/turin_n613_xeno_qn_genefilter_opt2_KMcalls_13052021.txt', sep="\t", header=T, stringsAsFactors = FALSE)

km <- km[km$KM.Subtypes != "Mixed",]
km$model <- substr(km$Sample_ID, 0, 7)
km <- km[grepl("LMX", km$Sample_ID),]

#meta <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_ctx', sep="\t", header=T)
#check all basali= yes!
### irino
mi <- merge(km, iri, by="model", all.x=T)
mi$class <- ifelse(is.na(mi$irino), 'Undetermined', ifelse(mi$irino < 35, 'Responder', "Non-Responder"))

pld <- mi %>% 
  group_by(KM.Subtypes,class) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

ggplot(pld, aes(x = factor(KM.Subtypes), y = perc*100, fill = factor(class))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "KM", y = "percent", fill = "Irinotecan 3Wd%") +
  theme_bw(base_size = 14)+scale_fill_manual(values=c('red','darkblue','grey'))
# 
# pld <- mix %>% 
#   group_by(KM.prediction,class) %>% 
#   summarise(count = n()) %>% 
#   mutate(perc = count/sum(count))
# 
# ggplot(pld, aes(x = factor(KM.prediction), y = perc*100, fill = factor(class))) +
#   geom_bar(stat="identity", width = 0.7) +
#   labs(x = "KMpred", y = "percent", fill = "Irinotecan 3Wd%") +
#   theme_bw(base_size = 14)+scale_fill_manual(values=c('red','darkblue','grey'))

table(mix$KM.Subtypes)
table(mix$KM.prediction)
### cetuxi
mi <- merge(km, ctx, by="model", all.x=T)
mi$class <- ifelse(is.na(mi$ctx), 'Undetermined', ifelse(mi$ctx < 35, 'Responder', "Non-Responder"))

pld <- mi %>% 
  group_by(KM.Subtypes,class) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

ggplot(pld, aes(x = factor(KM.Subtypes), y = perc*100, fill = factor(class))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "KM", y = "percent", fill = "Cetuxi 3Wd%") +
  theme_bw(base_size = 14)+scale_fill_manual(values=c('red','darkblue','grey'))


table(mix$KM.Subtypes)
table(mix$KM.prediction)
## mek inibitori


############## subset w.r.t what's in CRIS ##############