library(DESeq2)
library(ggplot2)
library(gridExtra)

pca_data <- snakemake@input[["pca"]]
leuco_data <- snakemake@input[["leuco"]]
figure <- snakemake@output[["fig"]]

#load('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/pc_manual_plots.Rdata')
load(pca_data)
what <- 'type'
p0 <- pc_cool(object = vsd, intgroup = what, pc1 = number1, pc2 = number2, plotfile=pca, legend=TRUE)
p1 <- pc_cool(object = vsd, intgroup = what, pc1 = 3, pc2 = 4, plotfile="meh.pdf", legend=FALSE)
what <- 'batch'
p2 <- pc_cool(object = vsd, intgroup = what, pc1 = 1, pc2 = 2, plotfile="meh.pdf", legend=TRUE)
p3 <- pc_cool(object = vsd, intgroup = what, pc1 = 3, pc2 = 4, plotfile="meh.pdf", legend=FALSE)
p4 <- pc_cool(object = vsd, intgroup = what, pc1 = 5, pc2 = 6, plotfile="meh.pdf", legend=FALSE)
p5 <- pc_cool(object = vsd, intgroup = what, pc1 = 7, pc2 = 8, plotfile="meh.pdf", legend=FALSE)
p6 <- pc_cool(object = vsd, intgroup = what, pc1 = 9, pc2 = 10, plotfile="meh.pdf", legend=FALSE)

#load("/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/leuco_score.Rdata")
load(leuco_data)
if ( color == "Leucocyte" | color == "CAF" | color == "Endothelial") {
  p7 <- ggplot(df, aes_string(x="PC1", y=y, color=color, shape=shape)) +
    labs(title="Scatter plot PC") +
    scale_color_distiller(type="seq", palette="OrRd", direction=1) +
    scale_shape_manual(values=c(0,1,2,3,7,16,17,8)) +
    geom_point() +
    scale_y_continuous(breaks=seq(-50, 70, 10)) +
    geom_hline(yintercept=0, color="black") +
    geom_vline(xintercept=0, color="black") +
    theme_bw()+unmute_theme #+
  # theme(axis.line=element_line(size=2, colour = "black"))
} else {
  p7 <- ggplot(df, aes_string(x="PC1", y=y, color=color, shape=shape)) +
    labs(title="Scatter plot PC") +
    scale_color_brewer(type='qual', palette="Dark2") +
    scale_shape_manual(values=c(0,1,2,3,7,16,17,8)) +
    geom_point() + 
    scale_y_continuous(breaks=seq(-50, 70, 10)) +
    geom_hline(yintercept=0, color="black") +
    geom_vline(xintercept=0, color="black") +
    theme_bw()+unmute_theme #+
  # theme(axis.line=element_line(size=2, colour = "black"))
}
lay <- rbind(c(1, 2), c(3, 4), c(5,6), c(7,8))
m <- grid.arrange(p0, p1, p2, p3, p4, p5, p6, p7, layout_matrix=lay)
ggsave(filename=figure, plot=m, dpi=300, width=8.3, height=11.7, units="in")
#savehistory('testPC.Rhistory')
