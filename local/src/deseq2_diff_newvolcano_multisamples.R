#print(snakemake@params[["threads"]]) # rid to name
library(DESeq2)
library(ggplot2)
library(ggrepel)

#register(MulticoreParam(as.numeric(snakemake@params[["threads"]])))
threads <- as.numeric(snakemake@params[["threads"]])
parallel <- FALSE
if (threads > 1) {
  library("BiocParallel")
  register(MulticoreParam(threads))
  parallel <- TRUE
}
alpha <- as.numeric(snakemake@params[["alpha"]])
lfc <- as.numeric(snakemake@params[["lfc"]]) # used only for volcano plots, the tsv printed lists all non NA results!

#print(snakemake@input[[1]]) # .RData
#print(snakemake@params[["class"]]) # which columns need to be compared
#print(snakemake@params[["nom"]]) # this vs
#print(snakemake@params[["den"]]) # this one
con <- c(snakemake@params[["factor"]], snakemake@params[["nom"]], snakemake@params[["den"]])
volcano <- snakemake@output[["volcano"]]
tsv <- snakemake@output[["tsv"]]
#load overwrites our snakemake object thus we need to put aside our parameters before.

save.image(paste0(tsv, "_DESeq.Rdata"))

load(snakemake@input[[1]])

res <- results(dds, alpha=alpha, contrast=con, parallel=parallel)
plot_volcano <- function(resnona, alpha, lfc, outfile, title) {
  
  # --- LOGICA DATI (TUA ORIGINALE) ---
  for (i in rownames(resnona)) {
    if (resnona[i,"log2FoldChange"] > lfc && resnona[i,"padj"] < alpha) {
      resnona[i, "sign"] <- "up"
    } else if (resnona[i,"log2FoldChange"] < -lfc && resnona[i,"padj"] < alpha) {
      resnona[i, "sign"] <- "down"
    } else {
      resnona[i,"sign"] <- "no_diff"
    }
  }

  resnona$signtot <- ifelse(abs(resnona$log2FoldChange) > lfc & resnona$padj < alpha, "both", 
                            ifelse(abs(resnona$log2FoldChange) > lfc, "LFC",
                                   ifelse(resnona$padj < alpha, "padj", "NS")))
  
  # --- CALCOLO DEL LIMITE Y MASSIMO ---
  # Cap i valori di padj troppo piccoli per evitare infiniti
  resnona$padj_capped <- pmax(resnona$padj, 1e-300)
  
  # Calcolo -log10 sui valori cappati
  log10_padj <- -log10(resnona$padj_capped)
  
  # Trovo il valore massimo
  max_y_val <- max(log10_padj[is.finite(log10_padj)], na.rm = TRUE)
  
  # Definisco il nuovo "soffitto" del grafico con più spazio
  y_limit_upper <- max_y_val * 1.3 
  
  # --- PARTE GRAFICA ---
  p <- ggplot(resnona, aes(log2FoldChange, -log10(padj_capped))) +
    geom_point(aes(col = sign), size=0.7, alpha=0.6) + 
    theme_bw() +
    scale_color_manual(values = c("blue", "#999999", "red"), drop=FALSE) + 
    ggtitle(title) +
    
    # QUI LA SOLUZIONE:
    # Imposto i limiti HARD dell'asse Y. 
    # c(0, y_limit_upper) forza l'asse ad arrivare fino a lassù.
    # La riga nera superiore sarà disegnata a y_limit_upper.
    scale_y_continuous(limits = c(0, y_limit_upper)) +
    
    coord_cartesian(clip = "off")

  # --- ETICHETTE ---
  nsign <- nrow(resnona[resnona$signtot=="both",])
  
  # Aggiungo 'nudge_y' per dare una spintarella in su alle etichette
  # Aggiungo 'force' per aumentare la repulsione
  repel_params <- list(
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.5,
    force = 2,              # Spinge via le etichette con più forza
    nudge_y = max_y_val * 0.05 # Sposta l'ancoraggio leggermente in su
  )

  if (nsign > 20) {
    p_final <- p + do.call(geom_text_repel, c(list(data=resnona[1:10,], aes(label=rownames(resnona)[1:10])), repel_params))
  } else {
    p_final <- p + do.call(geom_text_repel, c(list(data=resnona[resnona$signtot=="both",], aes(label=rownames(resnona[resnona$signtot=="both",]))), repel_params))
  }
  
  ggsave(outfile, plot = p_final, width = 8, height = 10)
}
resnona <- res[!is.na(res$pvalue) & !is.na(res$padj),]
resnona_df <- as.data.frame(resnona[order(resnona$padj),])
title <- trimws(strsplit(elementMetadata(res)[2,2], ":")[[1]][2])
plot_volcano(resnona_df, alpha, lfc, volcano, title)

write.table(resnona_df, file=tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

save.image(paste0(tsv, "_DESeq.Rdata"))