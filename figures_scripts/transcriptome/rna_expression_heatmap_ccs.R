library(tidyverse)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(ComplexHeatmap)
library(openxlsx)
library(InteractiveComplexHeatmap)
BiocManager::install("biomaRt")

setwd("C:/Users/dawid/Desktop/transcriptome")
setwd("/Users/dawid/sciebo/drArbeit/diverse/data/transcriptome")
setwd("/Users/dawid/sciebo/drArbeit/diverse/data/transcriptome")
getwd()

cts <- read.table("TPM.csv", 
                  header = TRUE, sep = ",", row.names = 1)

dim(cts)
#View(cts)

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensg=row.names(cts)

#ensg=c("ENSG00000147889","ENSG00000160072", "ENSG00000234396", "ENSG00000225972", "ENSG00000224315", "ENSG00000198744", "ENSG00000279928", "ENSG00000228037")

annot = getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'transcript_length', 'transcript_is_canonical', 'gene_biotype'),
              filters = 'ensembl_gene_id', 
              values = ensg, 
              mart = ensembl,useCache = FALSE,
              uniqueRows = TRUE)



annot$hgnc_symbol <- ifelse(annot$hgnc_symbol=="", annot$ensembl_gene_id, annot$hgnc_symbol)
annot_backup=annot
annot=annot_backup
# only canonical genes
annot=drop_na(annot)
# only coding 
annot=annot[annot$gene_biotype == "protein_coding", ]

annot_safe=annot
annot$transcript_is_canonical=NULL

annot=annot[!duplicated(annot$ensembl_gene_id),]
row.names(annot)=annot$ensembl_gene_id

df=merge(cts, annot, by="row.names", all=TRUE)
# all annoted?? Checken
df=drop_na(df)

row.names(df)=df$Row.names
#row.names(df)=paste(df$hgnc_symbol, df$Row.names, sep= "_")
df$Row.names=NULL
df$hgnc_symbol=NULL
df$ensembl_gene_id=NULL
df$gene_biotype=NULL
df$transcript_length.y=NULL


df=as.data.frame(df)

write.csv(df, "raw_reads_coding_only_CCS.csv",  row.names=TRUE)

rownames(df) <- df$Row.names
write.xlsx(df, 'name-of-your-excel-file.xlsx')

l= colnames(df)

l

df2[,c(1,3,2,4)]

df["ENSG00000157404",]

##########################
##### TPM in jupyter #####
##########################

tpm=read.csv("TPM_gist.csv", row.names=1)

dim(tpm)
head(tpm)
### tpm unter 1 als null???
# tpm[tpm < 1] <- 0   

tpm <- tpm[ rowSums(tpm) > 1, ]

log_tpm=tpm

log_tpm <- log2(tpm + 1 )

tpm["ENSG00000147889",]

######## selected gens ######
genes = c("SDHA","SDHB","SDHC","SDHD")
genes = c("CDKN2A","CDKN2B")


annot_select = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), 
              filters = 'hgnc_symbol', 
              values = genes, 
              mart = ensembl,useCache = FALSE)

selected=annot_select$ensembl_gene_id
heatmat_select = log_tpm[selected,]
rownames(annot_select)=annot_select$ensembl_gene_id

heatmat_select = merge(heatmat_select, annot_select, by="row.names")
rownames(heatmat_select)= heatmat_select$hgnc_symbol
heatmat_select$Row.names = NULL
heatmat_select$hgnc_symbol = NULL
heatmat_select$ensembl_gene_id = NULL

heatmat_select=as.matrix(heatmat_select)


title = c("")

hmcol <- colorRampPalette(brewer.pal(3, "OrRd"))(10) ## hmcol <- heat.colors
#pdf("PCA_CCS.pdf")
heatmap.2(heatmat_select, scale="row", Rowv = TRUE,  
          dendrogram="both", trace="none", margin=c(5,5), cexRow=1, cexCol=0.8, keysize=1, col = hmcol,
          labRow = NULL, Colv=TRUE)
#dev.off()


hmcol=c("#Ffffff","#E34A33")


a=colnames(heatmat_select)
a
b=c("loss (haploid)","deletion","diploid","diploid","diploid","diploid","diploid","diploid","loss (haploid)","loss (haploid)",
    "deletion","diploid","loss (haploid)","diploid","loss (haploid)","deletion","deletion","loss (haploid)","loss (haploid)","diploid",
    "diploid","loss (haploid)","deletion","deletion","diploid")

df=data.frame(a,b)

ha_column = HeatmapAnnotation(CNVkit_CDKN2A_status=df$b)


#pdf("PCA_CCS.pdf")
Heatmap(heatmat_select, name = "log2(tpm)", na_col = "black",
        column_title = title ,top_annotation = ha_column )
#dev.off()


ht = Heatmap(heatmat_select, name = "log2(tpm)", na_col = "black",
             column_title = title ,top_annotation = ha_column )
ht = draw(ht) # not necessary, but recommended

htShiny(ht)





#write.csv(round(t(heatmat_select), digits = 1), file = "/Users/dawid/sciebo/drArbeit/diverse/data/transcriptome/tpm_coding.csv", row.names = TRUE)

