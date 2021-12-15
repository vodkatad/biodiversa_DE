library(DESeq2)
library(tidyverse)


load(snakemake@input[[1]])
pca <- snakemake@output[["pca_cool"]]
number1 <- as.numeric(snakemake@wildcards[["pc1"]])
number2 <- as.numeric(snakemake@wildcards[["pc2"]])
what <- snakemake@wildcards[["what"]]


pc_cool <- function (object, intgroup = "condition", pc1=1, pc2 =2, ntop = 500, returnData = FALSE, plotfile) {
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
      geom_point(size = 2, alpha = 0.3) + xlab(paste0(colnames(d)[1],": ", round(percentVar[pc1] * 100), "% variance")) + 
      ylab(paste0(colnames(d)[2],": ", round(percentVar[pc2] * 100), "% variance"))+ coord_fixed() + theme_bw()
    
    ggsave(filename = plotfile, plot = p1) 
    #print(p1)
  }
vsd <- vst(dds, blind = FALSE)

pc_cool(object = vsd, intgroup = "type", pc1 = number1, pc2 = number2, plotfile=pca)
