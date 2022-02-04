array<- read.table('nc_cris_uarray_0.2.tsv' ,sep="\t", header=T, stringsAsFactors = F)
us <- read.table("cris_fpkm_lmx_nc_arm.tsv",sep="\t", header=T, stringsAsFactors = F)
us$cris <- gsub("-", "", us$cris, fixed=T)
m<- merge(array, us, by= "genealogy")

table(m$cris.x, m$cris.y)
table(us$cris)
table(array$cris)