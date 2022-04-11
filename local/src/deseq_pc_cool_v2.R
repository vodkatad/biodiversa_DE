library(DESeq2)
#library(tidyverse)
library(ggplot2)

load(snakemake@input[[1]])
pca <- snakemake@output[["pca_cool"]]
number1 <- as.numeric(snakemake@wildcards[["pc1"]])
number2 <- as.numeric(snakemake@wildcards[["pc2"]])
what <- snakemake@wildcards[["what"]]


pc_cool <- function (object, intgroup = "condition", pc1=1, pc2 =2, ntop = 500, returnData = FALSE, plotfile, legend=TRUE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    if (intgroup == "type") {
      group <- as.character(group)
      group <- as.factor(unlist(lapply(strsplit(group, "_", fixed=TRUE), function(x) {x[[1]][1]})))
    }
    d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group, 
                    intgroup.df, name = colnames(object))
    colnames(d)[1] <- paste0("PC", pc1)
    colnames(d)[2] <- paste0("PC", pc2)
    
    if (returnData) {
      attr(d, "percentVar") <- percentVar[pc1:pc2]
      return(d)
    }
    # alpha and size
    
    p1 <- ggplot(data = d, aes_string(x = colnames(d)[1], y = colnames(d)[2], color = "group")) + 
      geom_point(size = 0.5, alpha = 0.3) + xlab(paste0(colnames(d)[1],": ", round(percentVar[pc1] * 100), "% variance")) + 
      ylab(paste0(colnames(d)[2],": ", round(percentVar[pc2] * 100), "% variance"))+ coord_fixed() + theme
    if (! legend) {
      p1 <- p1 + theme(legend.position="none")
    }
    ggsave(filename = plotfile, plot = p1, dpi=300, height=120, width=120, units='mm') 
    #print(p1)
    return(p1)
  }
vsd <- vst(dds, blind = FALSE)


textSize <- 5
largerSize <- 7

theme <- theme_bw() +
theme(
	text = element_text(size = textSize, family='sans'),
	axis.title = element_text(size = largerSize),
	axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
	axis.text.y = element_text(size = textSize, color="black"),
	plot.title = element_text(size = largerSize, hjust = 0.5),
	legend.title = element_text(size=largerSize),
    legend.text = element_text(size=textSize)
)


p <- pc_cool(object = vsd, intgroup = what, pc1 = number1, pc2 = number2, plotfile=pca, legend=TRUE)

save.image('pc_manual_plots.Rdata')