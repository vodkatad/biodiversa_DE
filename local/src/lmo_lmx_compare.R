## Script to compare expression levels of paired models between LMO-LMX.
library(ggplot2)
library(dplyr)
# We will get two input files with expression levels
lmo_expr_f <- snakemake@input[["lmo"]]
lmx_expr_f <- snakemake@input[["lmx"]]

#### wip ###
#setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')
#lmo_expr_f <- 'LMO_BASALE-CYP2W1_tmm.tsv'
#slmx_expr_f <- 'LMX_BASALE-CYP2W1_tmm.tsv'
############

# We will produce two output files:
scatter_plot_f <- snakemake@output[["scatter"]]
common_quartiles_f <- snakemake@output[["quartiles"]]

# Load files
lmo_expr_df <- read.table(lmo_expr_f, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")
lmx_expr_df <- read.table(lmx_expr_f, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")
colnames(lmo_expr_df) <- c('genealogy','expr')
colnames(lmx_expr_df) <- c('genealogy','expr')

# Log2(value+1)

# R "base" version
# lmo_expr_df$log2_expr <- log2(lmo_expr_df$expr+1)
# lmx_expr_df$log2_expr <- log2(lmx_expr_df$expr+1)
lmo_expr_df <- mutate(lmo_expr_df, log2_expr = log2(expr+1))
lmx_expr_df <- mutate(lmx_expr_df, log2_expr = log2(expr+1))

# Compute averages at the model level
lmo_expr_df <- mutate(lmo_expr_df, model = substr(genealogy, 0, 7))
lmx_expr_df <- mutate(lmx_expr_df, model = substr(genealogy, 0, 7))

# Merge the two data.frames and do the scatter plot + correlation
average_model <- function(model, df, col_value, col_average) {
  my_df <- df[df[,col_value] == model,]
  return(mean(my_df[, col_average]))
}

# models <- unique(lmo_expr_df$model)
# averaged_lmo <- data.frame(row.names=models, average=rep(NA, length(models)))
# for (i in seq(1, length(models))) {
#   averaged_lmo[i,'average'] <- average_model(models[i], lmo_expr_df, 'model','log2_expr')
# }

averaged_lmo <- as.data.frame(sapply(unique(lmo_expr_df$model), average_model, lmo_expr_df, 'model','log2_expr'))
colnames(averaged_lmo) <- c('log2_expr')
averaged_lmx <- as.data.frame(sapply(unique(lmx_expr_df$model), average_model, lmx_expr_df, 'model','log2_expr'))
colnames(averaged_lmx) <- c('log2_expr')

merged_expr <- merge(averaged_lmo, averaged_lmx, by="row.names")
colnames(merged_expr) <- c('model','lmo','lmx')

pearson <- cor.test(merged_expr$lmo, merged_expr$lmx)
correlation <- pearson$estimate
pvalue <- pearson$p.value
ggplot(data = merged_expr, aes(x=lmx, y=lmo)) + geom_point() + theme_bw() +
labs(caption=paste0(round(correlation, digits=3), ', pval=', formatC(pvalue, format = "e", digits = 3))) +
theme(text=element_text(size=15))

ggsave(scatter_plot_f)

# Prepare the output file with n. of models and common quartiles
res <- data.frame(class="tot", num=nrow(merged_expr), stringsAsFactors = FALSE)

# Count the n. of common models between quartiles
n_quartiles <- ceiling(nrow(merged_expr) / 4)

lmo_ordered <- merged_expr[order(merged_expr$lmo),]
lmx_ordered <- merged_expr[order(merged_expr$lmx),]

last_index <- 1 # variable to track up to where we got quartiles
for (i in seq(1,4)) {
  next_last_index <- last_index+n_quartiles
  if (next_last_index > nrow(merged_expr)) {
    next_last_index <- nrow(merged_expr)
  }
  x <- lmx_ordered[seq(last_index, next_last_index),]
  o <- lmo_ordered[seq(last_index, next_last_index),]
  common <- length(intersect(x$model, o$model))
  res[i+1,'class'] <- paste0('q', i)
  res[i+1,'num'] <- common
  last_index <- next_last_index+1 # Update the starting point for next interation
}

# Output files
write.table(res, file=common_quartiles_f, sep="\t", col.names = TRUE, row.names= FALSE, quote = FALSE)

