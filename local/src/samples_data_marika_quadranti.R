## samples data marika for global deg
library(tidyverse)

c <- snakemake@input[["counts"]]
kras <- snakemake@input[["kras_ref"]]
of <- snakemake@output[["o"]]
tipo <- snakemake@wildcards[["type"]]

#c <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Marika_RNA_26/counts_ok.gz"
c <- read.table(c, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
c$Geneid <- NULL

samples <- as.data.frame(matrix(ncol = 1, nrow = length(colnames(c))))
samples$V1 <- colnames(c)
samples$smodel <- substr(samples$V1, 1, 7)
samples$treat <- str_extract(samples$V1, "(?<=_)[^_]+(?=_)")
samples$replicate <- str_extract(samples$V1, "[^_]+$")
rownames(samples) <- samples$V1
samples$V1 <- NULL
samples$replicate <- paste0("R", samples$replicate)

#d <- read.table('/mnt/cold1/snaketree/prj/pdxopedia/local/share/data/kras_dep2026.tsv', sep="\t", header=T)
ggplot(d, aes(x = rasdep, y=tram, label=smodel))+geom_point()+
  geom_text(vjust = -1, size = 2)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.25, linetype = "dashed")

d <- read.table(kras, sep="\t", header=T)
sens_group <- d[d$rasdep < 0.25 & d$tram < 0.25,]
other_group <- d[d$rasdep > 0.25 & d$tram > 0.25,]

samples <- samples %>%
  mutate(group = case_when(
    smodel %in% sens_group$smodel ~ "sens",
    smodel %in% other_group$smodel ~ "other",
    TRUE ~ NA_character_
  ))

samples <- samples[!is.na(samples$group), ]

samples <- samples %>% filter(group == tipo)

write.table(samples, file=of, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)