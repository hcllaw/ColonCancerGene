setwd("~/Desktop/ColonCancerGene")
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
gse33113 <- getGEO("GSE33113")[[1]]
healthy= c('GSM1100477','GSM1100478','GSM1100479','GSM1100480','GSM1100481','GSM1100482')
gse33113 <- gse33113[,!(sampleNames(gse33113) %in% healthy)]
x <- exprs(gse33113)
dim(x)
#Dist_fgene = dist(full_gene_data)
#cluster_gene = hclust(Dist_fgene, method = "average")
source("https://bioconductor.org/biocLite.R")
biocLite("frma")
biocLite("barcode")
library(frma)
library(barcode)
frma(x)
object = barcode(x)
biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
x[1:5,1:4]
x[1,]
pamoutput <- read.csv('desousaemelo_146genes_pam.csv',header=TRUE)
head(pamoutput)
nrow(pamoutput)
x146 <- x[as.character(pamoutput$identifier),]
