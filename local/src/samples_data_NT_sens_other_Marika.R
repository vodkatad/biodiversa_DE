## samples data for NT sens other
library(tidyverse)

c <- snakemake@input[["counts"]]
kras <- snakemake@input[["kras_ref"]]
of <- snakemake@output[["o"]]
tipo <- snakemake@wildcards[["type"]]

#c <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/Marika_global/counts_ok.tsv.gz"
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
samples <- samples %>% filter(treat == "NT")

#d <- read.table('/mnt/cold1/snaketree/prj/pdxopedia/local/share/data/kras_dep2026.tsv', sep="\t", header=T)
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
samples$treat <- NULL
write.table(samples, file=of, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)