### load multipe tsv
# i created a directory with all copied files for differential expressed genes with
# mkdir diff_genes; find . -name "*deseq2.tsv" -exec cp -r {} diff_genes \;

##Read files named
filenames <- list.files(path="/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/diff_genes/",
                        pattern="*.tsv")

##Create list of data frame names without the extra part 
names <-substr(filenames,1,7)

###Load all files
for(i in names){
  filepath <- file.path("/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2_2nd/diff_genes/",paste(i,"_geno_cutoff0.05-NE.vs.N2.deseq2.tsv",sep=""))
  assign(i, read.table(filepath, quote="", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
}  

#prova <- CRC0059
# for (i in rownames(prova)) {
#   if (prova[i, "log2FoldChange"] >0){
#     prova[i, "DEG"] <- "UP"
#   } else {
#     prova[i, "DEG"] <- "DOWN"
#   }
# }
# prova_crc <- table(prova$DEG)

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

numberdiffgenes <- function(df) {
  df <- df %>% filter(padj<0.05)
  for (i in rownames(df)) {
    if (df[i, "log2FoldChange"] >0){
      df[i, "DEG"] <- "UP"
    } else {
      df[i, "DEG"] <- "DOWN"
    }
  }
  t <- table(df$DEG)
  return(t)
}

CRC0059 <- numberdiffgenes(CRC0059)
CRC0065 <- numberdiffgenes(CRC0065)
CRC0078<-numberdiffgenes(CRC0078)
CRC0080<-numberdiffgenes(CRC0080)
CRC0096<-numberdiffgenes(CRC0096)
CRC0116<-numberdiffgenes(CRC0116)
CRC0148<-numberdiffgenes(CRC0148)
CRC0152<-numberdiffgenes(CRC0152)
CRC0161<-numberdiffgenes(CRC0161)
CRC0196<-numberdiffgenes(CRC0196)
CRC0277<-numberdiffgenes(CRC0277)
CRC0291<-numberdiffgenes(CRC0291)
CRC0316<-numberdiffgenes(CRC0316)
CRC0322<-numberdiffgenes(CRC0322)
CRC0327<-numberdiffgenes(CRC0327)
CRC0399<-numberdiffgenes(CRC0399)
CRC0456<-numberdiffgenes(CRC0456)
CRC0464<-numberdiffgenes(CRC0464)
CRC0515<-numberdiffgenes(CRC0515)
CRC0534<-numberdiffgenes(CRC0534)
CRC0542<-numberdiffgenes(CRC0542)
CRC1239<-numberdiffgenes(CRC1239)
CRC1272<-numberdiffgenes(CRC1272)
CRC1278<-numberdiffgenes(CRC1278)
CRC1307<-numberdiffgenes(CRC1307)
CRC1314<-numberdiffgenes(CRC1314)
CRC1331<-numberdiffgenes(CRC1331)
CRC1430<-numberdiffgenes(CRC1430)
CRC1502<-numberdiffgenes(CRC1502)
CRC1588<-numberdiffgenes(CRC1588)
CRC1729<-numberdiffgenes(CRC1729)

# tables <- as.data.frame(matrix(ncol = 3))
# colnames(tables) <- c("CRC", "DOWN", "UP")
# for(i in dfs) {
#   for (j in rownames(i)) {
#     if (i[j, "log2FoldChange"] >0){
#       i[j, "DEG"] <- "UP"
#     } else {
#       i[j, "DEG"] <- "DOWN"
#     }
#   }
#   
#   tables[,"DOWN"] <- rbind(table(i$DEG)[1])
#   tables[,"UP"] <- rbind(table(i$DEG)[2])
#   tables[,"CRC"] <- ""
# }

tables <- as.data.frame(matrix(ncol = 3, nrow = 31))
colnames(tables) <- c("CRC", "DOWN", "UP")
tables$CRC <- names
#rownames(tables) <- tables$CRC

lista_tabelle <- list(CRC0059,CRC0065,CRC0078,CRC0080,CRC0096,CRC0116,CRC0148,CRC0152,CRC0161,
CRC0196,CRC0277,CRC0291,CRC0316,CRC0322,CRC0327,CRC0399,CRC0456,CRC0464,
CRC0515,CRC0534,CRC0542,CRC1239,CRC1272,CRC1278,CRC1307,CRC1314,CRC1331,
CRC1430,CRC1502,CRC1588,CRC1729)

tables <- as.data.frame(matrix(ncol = 3, nrow = 31))
colnames(tables) <- c("CRC", "DOWN", "UP")
tables$CRC <- names
rownames(tables) <- tables$CRC
tables$CRC <- NULL

tables["CRC0059",] <- lista_tabelle[[1]]
tables["CRC0065",] <- lista_tabelle[[2]]
tables["CRC0078",]<-lista_tabelle[[3]]
tables["CRC0080",]<-lista_tabelle[[4]]
tables["CRC0096",]<-lista_tabelle[[5]]
tables["CRC0116",]<-lista_tabelle[[6]]
tables["CRC0148",]<-lista_tabelle[[7]]
tables["CRC0152",]<-lista_tabelle[[8]]
tables["CRC0161",]<-lista_tabelle[[9]]
tables["CRC0196",]<-lista_tabelle[[10]]
tables["CRC0277",]<-lista_tabelle[[11]]
tables["CRC0291",]<-lista_tabelle[[12]]
tables["CRC0316",]<-lista_tabelle[[13]]
tables["CRC0322",]<-lista_tabelle[[14]]
tables["CRC0327",]<-lista_tabelle[[15]]
tables["CRC0399",]<-lista_tabelle[[16]]
tables["CRC0456",]<-lista_tabelle[[17]]
tables["CRC0464",]<-lista_tabelle[[18]]
tables["CRC0515",]<-lista_tabelle[[19]]
tables["CRC0534",]<-lista_tabelle[[20]]
tables["CRC0542",]<-lista_tabelle[[21]]
tables["CRC1239",]<-lista_tabelle[[22]]
tables["CRC1272",]<-lista_tabelle[[23]]
tables["CRC1278",]<-lista_tabelle[[24]]
tables["CRC1307",]<-lista_tabelle[[25]]
tables["CRC1314",]<-lista_tabelle[[26]]
tables["CRC1331",]<-lista_tabelle[[27]]
tables["CRC1430",]<-lista_tabelle[[28]]
tables["CRC1502",]<-lista_tabelle[[29]]
tables["CRC1588",]<-lista_tabelle[[30]]
tables["CRC1729",]<-lista_tabelle[[31]]

tables$CRC <- rownames(tables)

order_dep <- read_xlsx("/mnt/cold1/snaketree/prj/DE_RNASeq/local/share/data/tcf7l2_order_def.xlsx")
order_dep <- order_dep[-13,]

tables_order <- tables[order(match(tables$CRC, order_dep$case)), ]

melt_data <- melt(tables_order, id = c("CRC"))
melt_data$CRC <- as.factor(melt_data$CRC)

ggplot(melt_data, aes(x = fct_inorder(CRC), y = value, fill = variable, colour = variable)) + 
  geom_bar(stat = "identity", position = "dodge")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
