library(ggplot2)
library(reshape)
cris_candiolo <- '~/cris_tmm_0.2_classes_lmx_basali_models_ns.tsv'
cris_charles <- '/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/CharlesRiver/cris_tmm_0.2_classes.tsv'

#cms_candiolo TODO
cms_charles <- '/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/CharlesRiver/vsd_CMScaller.tsv'
cms_caller_f <- as.data.frame(t(data.frame(CMS1=77, CMS2=136, CMS3=76, CMS4=79, NC=183)))
colnames(cms_caller_f)[1] <- 'Freq'

cms_charles_df <- read.table(cms_charles, sep="\t", header=TRUE, stringsAsFactors = F)
cms_charles_df[is.na(cms_charles_df$prediction), 'prediction'] <- "NC"
cms_charles_f <- as.data.frame(table(cms_charles_df$prediction)) 
#cms_charles_f$Freq <- cms_charles_f$Freq/sum(cms_charles_f$Freq)
cms_charles_f$x <- 'CharlesRiver'

#cms_caller_f$Freq <- cms_caller_f$Freq/sum(cms_caller_f$Freq)
cms_caller_f$x <- 'CMScaller'
cms_caller_f <- cbind(cms_charles_f$Var1, cms_caller_f)
colnames(cms_caller_f)[1] <- 'Var1'

df <- rbind(cms_caller_f, cms_charles_f)
colnames(df)[1] <- 'CMS'
colnames(df)[2] <- 'N'
#ggplot(data=df, aes(x=x, y=N, fill=CMS))+geom_col(position="stack")+theme_bw()+theme(text = element_text(size = 20))


dft <- cast(data=df, "CMS~x", value="N")
rownames(dft) <- dft$CMS
dft$CMS <- NULL
chisq.test(as.matrix(as.data.frame(dft)))


###############
cris_candiolo_df <- read.table(cris_candiolo, sep="\t", header=TRUE)
cris_charles_df <- read.table(cris_charles, sep="\t", header=TRUE)

df_candiolo <- as.data.frame(table(cris_candiolo_df$cris)) # TODO ADD NA to cris?
df_candiolo$x <- 'Candiolo'
#df_candiolo$Freq <- df_candiolo$Freq/sum(df_candiolo$Freq)
df_charles <- as.data.frame(table(cris_charles_df$cris))
#df_charles$Freq <- df_charles$Freq/sum(df_charles$Freq)
df_charles$x <- 'CharlesRiver'
df <- rbind(df_candiolo, df_charles)
colnames(df)[1] <- 'CRIS'
colnames(df)[2] <- 'N'
#ggplot(data=df, aes(x=x, y=N, fill=CRIS))+geom_col(position="stack")+theme_bw()+theme(text = element_text(size = 20))


dft <- cast(data=df, "CRIS~x", value="N")
dft <- dft[dft$CRIS!="NS", ] 
rownames(dft) <- dft$CRIS
dft$CRIS <- NULL
chisq.test(as.matrix(as.data.frame(dft)))

