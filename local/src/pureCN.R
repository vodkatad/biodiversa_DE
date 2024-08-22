## purcone per gli amici (purecn)

library("PureCN")

reference.file <- "/mnt/trcanmed/snaketree/task/annotations/dataset/gnomad/GRCh38.d1.vd1.fa"
bed.file <- "/mnt/trcanmed/snaketree/prj/snakegatk/local/share/data/xgen-exome-hyb-panel-v2-probes-hg38.bed"
mappability.file <- "/mnt/trcanmed/snaketree/prj/snakegatk/local/share/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
intervals <- import(bed.file)
mappability <- import(mappability.file)
setwd("/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/")
preprocessIntervals(intervals, reference.file, mappability = mappability, output.file = "intervalli_done")

bam.file <- "/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/align/markedDup_CRC1599NLH0000000000D03000V2.sorted.bam"
bam.file.tumor <- "/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/align/markedDup_CRC1599PRX0A02002TUMD03000V2.sorted.bam"
intervalli_done <- "/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/intervalli_done"
cores <- 6
bamcoverage_normal <- calculateBamCoverageByInterval(bam.file = bam.file, interval.file = intervalli_done, output.file = "bam_coverage_normal")
bamcoverage_tumor <- calculateBamCoverageByInterval(bam.file = bam.file.tumor, interval.file = intervalli_done, output.file = "bam_coverage_tumor")
