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
X146data = as.data.frame(t(x146))
gene_model = Mclust(X146data, G = 3)
plot(gene_model, X146data,what = "BIC")
coordProj(X146data, dimens = c(5,7), what ="classification", classification = gene_model$classification, parameters = gene_model$parameters)
dim(x146)
gene_name = x146[,1]
library(gtools)
ylabels = as.vector(gene_model$classification)
library(pamr)
gene146data = list(x=as.matrix(as.data.frame(x146)),y=as.matrix(as.data.frame(ylabels)))
gene.train = pamr.train(gene146data)
gene.results = pamr.cv(gene.train, gene146data)
pamr.plotcv(gene.results)
pamr.confusion(gene.results, threshold=0.05)
pamr.plotcen(gene.train, gene146data, threshold=0.1)
par(mar=c(1,1,1,1))
pamr.geneplot(gene.train, gene146data, threshold=0.1)
pamr.listgenes(gene.train, gene146data, threshold=0.1)
khan.scales <- pamr.adaptthresh(gene.train)




library(mclust)
gene_class = MclustDA( as.data.frame(t(x146[,1:60])), factor(as.vector((ylabels[1:60]))))




