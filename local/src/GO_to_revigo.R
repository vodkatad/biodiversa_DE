library(tidyverse)
library(ggplot2)

GO <- snakemake@input[["GO_r"]]
rev <- snakemake@output[["revigo"]]

go <- read.table(GO, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go <- go %>% filter(p.adjust < 0.05)
go <- go[order(go$p.adjust),]
go <- go[,c("ID","p.adjust")]

write.table(go, rev, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)