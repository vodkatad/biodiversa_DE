odx <- read.table(gzfile('/mnt/cold1/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_t_mean_gene_genealogyall.tsv.gz'), sep="\t", header=T)
odo <- read.table(gzfile('/mnt/cold1/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_t_mean_gene_genealogyall.tsv.gz'), sep="\t", header=T)

w <- intersect(colnames(dx), colnames(do))

dx <- odx[, w]
do <- odo[, w]

dxave <- rowMeans(dx)
#keep <- dxave[dxave > quantile(dxave)[4]]
keep <- dxave[dxave > median(dxave)]

sds <- apply(dx, 1, sd)
osds <- sds[order(sds)]

keepn <- intersect(names(keep), names(tail(osds, n=3000)))

keepn2 <- names(tail(osds, n=1000))

any('H_CES1'  %in% keepn2)

dx <- dx[keepn,]
do <- do[keepn,]

cc <- cor(t(dx), t(do), method="spearman") # TODO add method="spearman"

cors <- diag(cc)

cis <- as.numeric(unlist(cc[cors]))
summary(cis)
cc[rownames(cc)=="H_CES1", colnames(cc)=="H_CES1"]
hist(cors)

#> summary(cis)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.650127 -0.070789  0.001202  0.001802  0.073573  0.885784 

plot(unlist(dx['H_OLFM4',]),unlist(do['H_OLFM4',]))


#> cc[rownames(cc)=="H_CES1", colnames(cc)=="H_CES1"]
#[1] 0.7376242