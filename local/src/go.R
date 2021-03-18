library("org.Hs.eg.db",  quietly=TRUE)
library(topGO,  quietly=TRUE)
# A two column file with a gene id and a class
classesf <- snakemake@input[["classes"]]
# A single column with the universe
universef <- snakemake@input[["universe"]]
outfile <- snakemake@output[[1]]
id <- snakemake@params[["ids"]]
ontologies <- snakemake@params[["onto"]]
debug <- snakemake@params[["debug"]]

classes <- read.table(classesf, header=FALSE, sep="\t")
universe <- read.table(universef, header=FALSE, sep="\t")

goenrich <- function(ontology, namedgenes, id) {
    GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = namedgenes,
              geneSel = function(x) {x!=0},
              nodeSize = 10,
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = id)
    allGO = usedGO(object = GOdata) 
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(allGO), numChar=1000)
    allRes$pnom <- score(resultFisher)
    #- as.numeric(allRes$classicFisher)
    allRes$classicFisher <- NULL
    allRes
}

goenrichall <- function(genes, universe, ontologies, id) {
    namedgenes <- rep(0, nrow(universe))
    names(namedgenes) <- universe[,1]
    namedgenes[names(namedgenes) %in% genes] <- 1
    allont <- lapply(ontologies, goenrich, namedgenes, id)
    garbage <- lapply(seq(1,length(ontologies)), function(i) {allont[[i]]$ont <<- ontologies[[i]]})
    allontdf <- do.call(rbind, allont)
    allontdf
}

all_classes <- unique(classes[,1])
all <- lapply(all_classes, function(x) { goenrichall(classes[classes[,1]==x,2], universe, ontologies, id)  })
garbage <- lapply(seq(1,length(all_classes)), function(i) {all[[i]]$class <<- all_classes[[i]]})
save.image('pippo.Rdata')
resdf <- do.call(rbind, all)
resdf$padj <- p.adjust(resdf$pnom, method="BH")
write.table(resdf, gzfile(outfile), quote=FALSE, row.names =  FALSE, col.names = TRUE, sep="\t")


if (debug == "yes") {
  save.image(file=paste0(outfile,'.debug','.RData'))
}
