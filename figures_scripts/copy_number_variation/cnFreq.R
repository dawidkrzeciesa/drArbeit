library("GenVisR")
packageVersion("GenVisR")
library(ggplot2)

#setwd("/home/rtx/Schreibtisch/test_cn")



cnData=read.table(file = "/Volumes/Dokumente/drArbeit_backup/fertige_analyzen/cnv/CCS_cnv_call_tumor_content/merged.cns.tsv", header=TRUE)
# cnFreq(cnData, genome="hg38", CN_low_cutoff = 1.8, CN_high_cutoff = 2.1,facet_lab_size = 5, out = "plot")


layer1 <- theme(panel.spacing = unit(0, "lines"))
cnFreq(cnData, genome="hg38", CN_low_cutoff = 1.8, CN_high_cutoff = 2.1,facet_lab_size = 5, plotLayer=layer1)


df=cnFreq(cnData, genome="hg38", CN_low_cutoff = 1.8, CN_high_cutoff = 2.1,facet_lab_size = 5, out = "data")
df2=do.call(rbind.data.frame, df)


sapply(cnData, typeof)

# highlight cdkn2a locus

layer1 <- geom_vline(xintercept=c(21981526))
cnFreq(cnData, genome="hg38", plotChr="chr9", plotLayer=layer1)

cnFreq(cnData, genome="hg38", plotChr="chr9")
