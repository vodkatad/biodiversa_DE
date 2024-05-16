library(pheatmap)
set.seed(42)
nrandom <- rnorm(150, mean=0, sd=2)
decenter <- nrandom[nrandom>-0.5]
decenter <- c(decenter, rep(0, 8))
d <- matrix(decenter, ncol=5)

minv <- -4
maxv <- 4
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", "#FFFFFF",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette)

#pheatmap(d, cluster_rows = F, cluster_cols=F,
#        breaks = seq(-4, 4, by=1))