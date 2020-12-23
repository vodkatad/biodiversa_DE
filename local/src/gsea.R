library(WriteXLS)
library(fgsea)
library(ggplot2)

input <- snakemake@input[["tsv"]]
pathways_data <- snakemake@input[["pathways"]]
outdir <- snakemake@output[["outdir"]]
outtable <- snakemake@output[["outtable"]]
outtableall <- snakemake@output[["outtableall"]]
cores <- snakemake@params[["cores"]]
save.image(snakemake@params[["save"]])
load(pathways_data)
run_all_gsea <- function(de_data, pathways) {
  rnked_de <- de_data[order(de_data$sort),]
  rnk <- rnked_de$sort
  names(rnk) <- rnked_de$geneid
  res <- lapply(seq(1, length(pathways)), run_gsea, rnk, unique(de_data$name))
  do.call("rbind", res)
}

run_gsea <- function(pathway_i, rnk, name) {
  df <- fgsea(pathways[[pathway_i]], rnk, minSize=15, maxSize=1000, nperm=10000, nproc=cores) # 10000
  df <- as.data.frame(df)
  if (nrow(df) > 0) {
    df$gs <- names(pathways)[[pathway_i]]
    df$name <- name
  }
  df$leadingEdge <- NULL
  df
}

# already in the loaded rdata:
#pathways <- list(Mm.H, Mm.c2, Mm.c3, Mm.c4, Mm.c6, Mm.c7)
#names(pathways) <- c("hallmark", "curated", "motif", "computational", "oncogenic", "immunologic")
#pathways <- pathways[c(1)] # we keep only hallmark
all <- read.table(input, header=TRUE, sep="\t")

dir.create(outdir)
setwd(outdir)
#allres <- by(all, all$name, run_all_gsea, pathways)
# inefficient but simpler
allres <- lapply(levels(all$name), function(x) {vs <-all[all$name==x,]; run_all_gsea(vs, pathways) })
res <- do.call("rbind", allres)
res$padj_multi <- p.adjust(res$pval, method="BH")
sign <- res[res$padj_multi < 0.05,]
write_xls <- function(data) {
  name <- unique(data$name)
  data$padj <- NULL
  data$nMoreExtreme <- NULL
  data$name <- NULL
  data$size <- NULL
  data <- data[order(data$pval),]
  WriteXLS(data, ExcelFileName=paste0(paste(name, sep="_"), ".gsea.xls"), row.names=FALSE, col.names=TRUE)
}

garbage <- by(res, res$name, write_xls)

plot_es <- function(row, all_de_data) {
  #name <- as.character(unique(all_de_data$name))
  #name2 <- levels(sign$name)[as.numeric(row[9])]
  #cat(paste0(row[8], row[1], row[9], ".pdf", collapse="_"))
  de_data <- all_de_data[as.character(all_de_data$name) == row[9],]
  rnked_de <- de_data[order(de_data$sort),]
  rnk <- rnked_de$sort
  names(rnk) <- rnked_de$geneid
  #cat(head(pathways[[as.character(row[8])]][[as.character(row[1])]]))
  #cat("uppa\nappa")
  #cat(head(names(rnk)))
  plotEnrichment(pathways[[as.character(row[8])]][[as.character(row[1])]],rnk) +labs(title=row[1])
  ggsave(paste0(c(row[8], row[1], row[9], "pdf"), collapse="."))
}

sign <- sign[order(sign$name, sign$pval),]
#to_plot <- lapply(levels(sign$name), function(x) { w <- sign[sign$name==x,];  head(w, n=20)})
#top <- do.call("rbind", to_plot)
#garbage <- apply(top, 1, plot_es, all)
#save.image(snakemake@params[["save"]])
setwd("../")
table <- data.frame(path=sign$pathway, NES=sign$NES, vs=sign$name, pval=sign$pval, padj=sign$padj_multi)
write.table(table, outtable, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(res, outtableall, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
