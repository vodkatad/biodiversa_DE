library(tidyverse)
library(readxl)

samples_f <- snakemake@input[["samples_or"]] 
order_dep_f <- snakemake@input[["order_or"]] 
meta <- snakemake@output[["meta"]]
type <- snakemake@wildcards[['geno']]

#samples_f <- read_xlsx("/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/TCF7L2_2nd/samples_names.xlsx")
samples <- read_xlsx(samples_f)
samples$...4 <- NULL
samples$...5 <- NULL
samples <- samples[-c(133:171),]
samples$id <- samples$REPLICATES
samples$model <- substr(samples$REPLICATES, 1,7)
samples$geno <- substr(samples$REPLICATES, 9, 10)
samples$replicates <- substr(samples$REPLICATES, 12, 13)
samples$sample <- NULL
samples$`tube name` <- NULL
samples$REPLICATES <- NULL

#order_dep_f <- "/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx"
order_dep <- read_xlsx(order_dep_f)
order_dep$quartile <- ntile(order_dep$co_comp_N2, 4)
names(order_dep)[names(order_dep)=="case"] <- "model"

samples_m <- merge(samples, order_dep, by="model")
samples_m <- samples_m[order(-samples_m$quartile),]
samples_m$co_comp_N2 <- NULL
rownames(samples_m) <- samples_m$id

for (i in rownames(samples_m)) {
  if (samples_m[i, "quartile"]==4){
    samples_m[i, "quartile"] <- "IND_4Q"
  } else if (samples_m[i, "quartile"]==1) {
    samples_m[i, "quartile"] <- "DEP_1Q"
  } else {
    samples_m[i, "quartile"] <- samples_m[i, "quartile"]
  }
}



samples_m <- samples_m %>% filter(geno == type)
samples_m$geno_resp <- paste0(samples_m$geno, "_", samples_m$quartile)
rownames(samples_m) <- NULL
samples_m <- samples_m[,c(2,1,3,4,5,6)]
samples_m$geno <- NULL
samples_m$quartile <- NULL

write.table(samples_m, file=meta, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)