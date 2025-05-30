library(citccmst)
summary(citccmst)
load(list.files(system.file("extdata", package = "citccmst"), full.names = TRUE))
## il file è già corretto per i symbols deprecati
## setdiff(symbol, vsd$geni) #MCUB #TMEM237
geni <- "/home/mferri/geni_marisa.xlsx"
geni <- read.xlsx(geni)
symbol <- geni$Gene.Symbol
geni <- geni[,c(2,3)]

citvalid.exp.annot <- data.frame(id=rownames(citvalid.exp.norm), stringsAsFactors=FALSE, row.names=rownames(citvalid.exp.norm))


vsd <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/vsd_fake.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(vsd) <- gsub("H_", "", rownames(vsd))
vsd$Gene.Symbol <- rownames(vsd)
vsd <- vsd %>% filter(Gene.Symbol %in% symbol)
vsd <- merge(vsd, geni, by="Gene.Symbol")
vsd$Gene.Symbol <- NULL
rownames(vsd) <- vsd$Probe.Set.ID
vsd$Probe.Set.ID <- NULL

citvalid.citccmst <- cit.assignCcmst(data = vsd,
                                     data.annot = citvalid.exp.annot,
                                     data.colId = "id",
                                     data.colMap = "id",
                                     citccmst.colMap = "Probe.Set.ID",
                                     dist.method = "dqda",
                                     plot = T)
table(citvalid.citccmst$citccmst.core)
table(citvalid.citccmst$citccmst.mixed)
samples <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/samples_data"
samples <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples$genealogy <- rownames(samples)

citvalid.citccmst$genealogy <- rownames(citvalid.citccmst)

res <- merge(samples, citvalid.citccmst, by="genealogy", all.y=TRUE)
res$sample <- NULL
res$batch <- NULL
res$type[is.na(res$type)] <- "no_deg"
ggplot(res, aes(x = citccmst.core, fill = type)) +
  geom_bar(position = "dodge")+scale_fill_manual(values = c("grey", "red", "blue"))

res <- res %>% filter(!type == "no_deg")
res <- res %>% filter(!citccmst.confidence == "OUTLIER")
ggplot(res, aes(x = citccmst.core, fill = type)) +
  geom_bar(position = "dodge")+scale_fill_manual(values = c("red", "blue"))
res$model <- substr(res$genealogy, 1, 7)

res2 <- res %>%
  group_by(model) %>%
  summarise(
    genealogy_diff = if_else(n_distinct(genealogy) > 1, "Different", "Same"),  
    citccmst.core = paste(unique(citccmst.core), collapse = ", ")
  )

res2 <- merge(res, res2, by="model")
res2 <- res2[!duplicated(res2$model),]
res <- res2
names(res)[names(res) == "citccmst.core.y"] <- "citccmst.core"
res$model <- NULL
res$citccmst.core.x <- NULL

df_aggregato <- res %>%
  group_by(citccmst.core, type) %>%
  summarise(count = n()) %>%  
  ungroup()

df_proporzioni <- df_aggregato %>%
  group_by(citccmst.core) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(df_proporzioni, aes(x = "", y = prop, fill = type)) +
  geom_bar(stat = "identity", width = 1) +  
  coord_polar("y") +  
  facet_wrap(~citccmst.core) + 
  theme_void() + 
  theme(legend.position = "bottom") 

#core_gruppo <- "C1"
#df <- res
rownames(res) <- res$genealogy
fisher_test_per_gruppo <- function(core_gruppo) {
  for (i in rownames(df)) {
    if (df[i, "citccmst.core"]==core_gruppo) {
      df[i, "gruppo_confronto"] <- core_gruppo
    } else {
      df[i, "gruppo_confronto"] <- "Altro"
    }
  }
  tabella <- table(df$gruppo_confronto, df$type)
  if (nrow(tabella) == 2 && ncol(tabella) == 2) {
    fisher_result <- fisher.test(tabella)
    return(data.frame(gruppo = core_gruppo, p_value = fisher_result$p.value, estimate = fisher_result$estimate))
  } else {
    return(data.frame(gruppo = core_gruppo, p_value = NA, estimate = NA))
  }
}

risultati_fisher <- do.call(rbind, lapply(unique(res$citccmst.core), fisher_test_per_gruppo))
rownames(risultati_fisher) <- risultati_fisher$gruppo

df_tot <- df_aggregato
df_tot <- cast(df_tot, type ~ citccmst.core)
rownames(df_tot) <- df_tot$type
df_tot$type <- NULL
df_tot <- as.data.frame(df_tot)
df_tot[is.na(df_tot)] <- 0
df_tot <- as.data.frame(t(df_tot))
df_tot <- as.matrix(df_tot)

fisher.test(df_tot)
df_tot <- as.data.frame(df_tot)
df_tot$gruppo <- rownames(df_tot)
df_tot <- merge(df_tot, risultati_fisher, by="gruppo")

#write.xlsx(df_tot, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/risultati_fisher_marisa.xlsx")
write.xlsx(df_tot, file="/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/risultati_fisher_marisa_collassomodello.xlsx")

#SD+PD -> PD
#perc < -50 is "PR"
#perc > 35 is "PD"
#middle SD

# rownames(samples) <- samples$sample_id_R
# for (i in rownames(samples)) {
#   if (samples[i,"X3WKS"] < -50) {
#     samples[i,"class"] <- "PR"
#   } else {
#     samples[i,"class"] <- "PD"
#   }
# }
# PD PR 
# 82  9 

## provo sui BASALI
# metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
# meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# meda_f$RNA_marker <- NULL
# meda_f$RNA_PC <- NULL
# meda_f$METHYL_L <- NULL
# meda_f$FRA_L <- NULL
# meda_f$w3_cetuxi <- NULL
# meda_f$w3_irino <- NULL
# lmx <- meda_f %>% filter(type %in% c("LMX_BASALE", "LMX_BASALE.1"))

# vsd_all <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
# vsd_all <- read.table(vsd_all, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# rownames(vsd_all) <- gsub("H_", "", rownames(vsd_all))
# vsd_all$Gene.Symbol <- rownames(vsd_all)
# vsd_all <- vsd_all %>% filter(Gene.Symbol %in% symbol)
# vsd_all <- merge(vsd_all, geni, by="Gene.Symbol")
# vsd_all$Gene.Symbol <- NULL
# rownames(vsd_all) <- vsd_all$Probe.Set.ID
# vsd_all$Probe.Set.ID <- NULL
# vsd_all <- as.data.frame(t(vsd_all))
# vsd_all$genealogy <- rownames(vsd_all)
# vsd_all <- vsd_all %>% filter(genealogy %in% lmx$sample_id_R)
# vsd_all$genealogy <- NULL
# vsd_all <- as.data.frame(t(vsd_all))
# 
# citvalid.citccmst_alllmx <- cit.assignCcmst(data = vsd_all,
#                                      data.annot = citvalid.exp.annot,
#                                      data.colId = "id",
#                                      data.colMap = "id",
#                                      citccmst.colMap = "Probe.Set.ID",
#                                      dist.method = "dqda",
#                                      plot = T)
# casi <- res
# casi <- casi[,c(1,2)]
# citvalid.citccmst_alllmx$genealogy <- rownames(citvalid.citccmst_alllmx)
# res_alllmx <- merge(citvalid.citccmst_alllmx, casi, by="genealogy")
# 
# ggplot(res_alllmx, aes(x = citccmst.core, fill = type)) +
#   geom_bar(position = "dodge")
