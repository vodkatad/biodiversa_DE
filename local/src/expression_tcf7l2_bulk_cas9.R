library(corrplot)

vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/TCF7L2/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)

s_n2 <- c("CRC0059", "CRC0322", "CRC0277", "CRC0515", "CRC1278", "CRC1502")
s_n2 <- paste0(s_n2, "_NKO_bulk")

s_cas9 <- c("CRC0059", "CRC0322", "CRC0277", "CRC0515", "CRC1278", "CRC1502")
s_cas9 <- paste0(s_cas9, "_CAS9")

s <- c(s_n2, s_cas9)

vsd <- vsd[,(colnames(vsd) %in% s)]

### filtering expression data: which sd genes but not clear outliers/not expressed genes
### filter not expressed genes
means_vsd <- apply(vsd, 1, mean)
med_vsd <- median(means_vsd)

### I obtain a data.table that contain only the genes  mean expression over the median
vsd_filtered <- as.data.frame(vsd[means_vsd > med_vsd,])

### I want to calculate the standard deviation
sds_vsd <- apply(vsd_filtered, 1, sd) 

### now we keep the top 10% variable genes 
sds_vsd <- sds_vsd[order(-sds_vsd)]
n_vsd <- length(sds_vsd)
keep_vsd <- head(sds_vsd, round(0.10*n_vsd))
keep_genes_vsd <- names(keep_vsd)
desd_vsd <- vsd_filtered[rownames(vsd_filtered) %in% keep_genes_vsd,]
res <- desd_vsd

c59 <- res[,grepl("CRC0059", colnames(res))]
c277 <- res[,grepl("CRC0277", colnames(res))]
c322 <- res[,grepl("CRC0322", colnames(res))]
c515 <- res[,grepl("CRC0515", colnames(res))]
c1278 <- res[,grepl("CRC1278", colnames(res))]
c1502 <- res[,grepl("CRC1502", colnames(res))]

c59$log2_59 <- log2((c59$CRC0059_NKO_bulk+1)/(c59$CRC0059_CAS9+1))
c277$log2_277 <- log2((c277$CRC0277_NKO_bulk+1)/(c277$CRC0277_CAS9+1))           
c322$log2_322 <- log2((c322$CRC0322_NKO_bulk+1)/(c322$CRC0322_CAS9+1))
c515$log2_515 <- log2((c515$CRC0515_NKO_bulk+1)/(c515$CRC0515_CAS9+1))
c1278$log2_1278 <- log2((c1278$CRC1278_NKO_bulk+1)/(c1278$CRC1278_CAS9+1))
c1502$log2_1502 <- log2((c1502$CRC1502_NKO_bulk+1)/(c1502$CRC1502_CAS9+1))

res_log2 <- as.data.frame(matrix(nrow = length(rownames(res)), ncol = 6))
rownames(res_log2) <- rownames(res)
colnames(res_log2) <- c("CRC0059", "CRC0277", "CRC0322", "CRC0515", "CRC1278", "CRC1502")

res_log2$CRC0059 <- c59$log2_59
res_log2$CRC0277 <- c277$log2_277
res_log2$CRC0322 <- c322$log2_322
res_log2$CRC0515 <- c515$log2_515
res_log2$CRC1278 <- c1278$log2_1278
res_log2$CRC1502 <- c1502$log2_1502

# cor59_277 <- cor.test(res_log2$CRC0059, res_log2$CRC0277)
# cor59_322 <- cor.test(res_log2$CRC0059, res_log2$CRC0322)
# cor59_515 <- cor.test(res_log2$CRC0059, res_log2$CRC0515)
# cor59_1278 <- cor.test(res_log2$CRC0059, res_log2$CRC1278)
# cor59_1502 <- cor.test(res_log2$CRC0059, res_log2$CRC1502)
# 
# cor277_322 <- cor.test(res_log2$CRC0277, res_log2$CRC0322)
# cor277_515 <- cor.test(res_log2$CRC0277, res_log2$CRC0515)
# cor277_1278 <- cor.test(res_log2$CRC0277, res_log2$CRC1278)
# cor277_1502 <- cor.test(res_log2$CRC0277, res_log2$CRC1502)
# 
# cor322_515 <- cor.test(res_log2$CRC0322, res_log2$CRC0515)
# cor322_1278 <-cor.test(res_log2$CRC0322, res_log2$CRC1278)
# cor322_1502 <- cor.test(res_log2$CRC0322, res_log2$CRC1502)
# 
# cor515_1278 <- cor.test(res_log2$CRC0515, res_log2$CRC1278)
# cor515_1502 <- cor.test(res_log2$CRC0515, res_log2$CRC1502)
# 
# cor1278_1502 <- cor.test(res_log2$CRC1278, res_log2$CRC1502)

Mprova = cor(res_log2)
prova <- corrplot(Mprova, method = 'number')

ggplot(res_log2, aes(x=CRC0277, y=CRC0515)) + 
  geom_point()+
  geom_smooth(method=lm)

ggplot(res_log2, aes(x=CRC1278, y = CRC1502))+geom_point()+geom_smooth(method = lm)
