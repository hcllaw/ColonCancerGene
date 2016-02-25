library(ConsensusClusterPlus)
library(mclust)

# set the working dir to the dir where you have gene146.csv file 
setwd("~/ColonCancerGene")

set.seed(174457795)

gene <- read.table("gene146.csv", header = T, sep = " ")
gene <- t(scale(t(gene)))
gene <- sweep(gene, 1, apply(gene,1,median,na.rm=T))
dim(gene)

# methods: kmeans, mclust, hierarchicalclust
methods <- c(1:3)

# ##################################################
# across clustering methods comparison
# ##################################################

clus_meth_cmp <- rep(numeric(0),3)
res1 <- ConsensusClusterPlus(gene, maxK=3, reps=500, pItem=1, pFeature=1,
                             clusterAlg="km",distance="pearson",seed=3,plot="png")
label1 <- res1[[3]][["consensusClass"]]
clustcard1 <- apply(matrix(1:3,1,3), 2, function(x) table(label1)[[x]])

res2 <- ConsensusClusterPlus(gene, maxK=3, reps=500, pItem=1, pFeature=1,
                             clusterAlg="hc",distance="pearson",seed=3,plot="png")
label2 <- res2[[3]][["consensusClass"]]
clustcard2 <- apply(matrix(1:3,1,3), 2, function(x) table(label2)[[x]])

res3 <- Mclust(t(gene), 3, "VVI")
label3 <- res3$classification
clustcard3 <- apply(matrix(1:3,1,3), 2, function(x) table(label3)[[x]])

d_mat <- matrix(0, nrow = dim(gene)[2], ncol = dim(gene)[2])
for (j in 1:dim(gene)[2])
  if (j !=dim(gene)[2])
    for (k in (j+1):dim(gene)[2])
      d_mat[j,k] <- as.integer(label1[j]==label1[k] && label2[j]!=label2[k])
D <- apply(matrix(1:3,1,3), 2, function(x) sum(d_mat[label1==x,label1==x]))
M <- apply(as.matrix(clustcard1), 1, function(x) as.integer(x*(x-1)/2))
clus_meth_cmp[1] <-  mean(D / M)


d_mat <- matrix(0, nrow = dim(gene)[2], ncol = dim(gene)[2])
for (j in 1:dim(gene)[2])
  if (j !=dim(gene)[2])
    for (k in (j+1):dim(gene)[2])
      d_mat[j,k] <- as.integer(label1[j]==label1[k] && label3[j]!=label3[k])
D <- apply(matrix(1:3,1,3), 2, function(x) sum(d_mat[label1==x,label1==x]))
M <- apply(as.matrix(clustcard1), 1, function(x) as.integer(x*(x-1)/2))
clus_meth_cmp[2] <- mean(D / M)


d_mat <- matrix(0, nrow = dim(gene)[2], ncol = dim(gene)[2])
for (j in 1:dim(gene)[2])
  if (j !=dim(gene)[2])
    for (k in (j+1):dim(gene)[2])
      d_mat[j,k] <- as.integer(label1[j]==label1[k] && label3[j]!=label3[k])
D <- apply(matrix(1:3,1,3), 2, function(x) sum(d_mat[label1==x,label1==x]))
M <- apply(as.matrix(clustcard1), 1, function(x) as.integer(x*(x-1)/2))
 mean(D / M)


d_mat <- matrix(0, nrow = dim(gene)[2], ncol = dim(gene)[2])
for (j in 1:dim(gene)[2])
  if (j !=dim(gene)[2])
    for (k in (j+1):dim(gene)[2])
      d_mat[j,k] <- as.integer(label3[j]==label3[k] && label2[j]!=label2[k])
D <- apply(matrix(1:3,1,3), 2, function(x) sum(d_mat[label3==x,label3==x]))
M <- apply(as.matrix(clustcard3), 1, function(x) as.integer(x*(x-1)/2))
clus_meth_cmp[3] <- mean(D / M)

# The vector clus_meth_cmp contains the values of table 1 in the report
clus_meth_cmp


# ##################################################
# within clustering methods comparison 
# ##################################################
sd_vals <- c(0, 0.05, 0.1, 0.5, 1, 1.5, 2)
methods_rate <- matrix(numeric(0),length(methods),length(sd_vals))

for (m in methods){

  if (m == 1){
    res <- ConsensusClusterPlus(gene, maxK=3, reps=500, pItem=1, pFeature=1,
                                clusterAlg="km",distance="pearson",seed=3,plot="png")
    label <- res[[3]][["consensusClass"]]
    clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(label)[[x]])
  }
  
  else if (m == 2){
    res <- ConsensusClusterPlus(gene, maxK=3, reps=500, pItem=1, pFeature=1,
                                clusterAlg="hc",distance="pearson",seed=3,plot="png")
    label <- res[[3]][["consensusClass"]]
    clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(label)[[x]])
  }
  
  else {
    res <- Mclust(t(gene), 3, "VVI")
    label <- res$classification
    clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(label)[[x]])
    length(label)
  }
  
  M <- apply(as.matrix(clustcard), 2, function(x) as.integer(x*(x-1)/2))
  noisy_gene <- list()
  noisy_res <- list()
  noisy_label <- list()
  
  for(s in 1:length(sd_vals)){
    rate <- c()
    m_rate <- c()
    for (i in 1:10){
      noisy_gene[[i]] <- t(scale(t(gene + rnorm(dim(gene)[1]*dim(gene)[2],0,sd_vals[s]))))
      noisy_gene[[i]] <- sweep(noisy_gene[[i]], 1, apply(noisy_gene[[i]],1,median,na.rm=T))
      d_mat <- matrix(0, nrow = dim(gene)[2], ncol = dim(gene)[2])

      if (m == 1){
        noisy_res[[i]] <- ConsensusClusterPlus(noisy_gene[[i]], maxK=3, reps=1000, pItem=1, pFeature=1,
                                               clusterAlg="km",distance="pearson",seed=3,plot="png")
        noisy_label <- noisy_res[[i]][[3]][["consensusClass"]]
        noisy_clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(noisy_label)[[x]])
      }

      else if (m == 2){
        noisy_res[[i]] <- ConsensusClusterPlus(noisy_gene[[i]], maxK=3, reps=1000, pItem=1, pFeature=1,
                                               clusterAlg="hc",distance="pearson",seed=3,plot="png")
        noisy_label <- noisy_res[[i]][[3]][["consensusClass"]]
        noisy_clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(noisy_label)[[x]])
      }

      else {
        noisy_res[[i]] <- Mclust(t(noisy_gene[[i]]), 3, "VVI")
        noisy_label <- noisy_res[[i]]$classification
        noisy_clustcard <- apply(matrix(1:3,1,3), 2, function(x) table(noisy_label)[[x]])
      }

      for (j in 1:dim(gene)[2])
        if (j !=dim(gene)[2])
          for (k in (j+1):dim(gene)[2])
            d_mat[j,k] <- as.integer(label[j]==label[k] && noisy_label[j]!=noisy_label[k])
      noisy_D <- apply(matrix(1:3,1,3), 2, function(x) sum(d_mat[label==x,label==x]))
      noisy_M <- apply(as.matrix(clustcard), 1, function(x) as.integer(x*(x-1)/2))
      rate <- c(rate, mean(noisy_D / noisy_M))
      cat("\n step", i, "\n")
    }
    methods_rate[m,s] <- mean(rate)
    cat("\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n", mean(rate),
        "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  } 
}


# ##################################################
# Code to get fig.2 in the report
# ##################################################
library(ggplot2)
library(reshape2)
result <- data.frame(sd_vals, t(methods_rate))
colnames(result) = c("sd","cons_kmeans","mclust","cons_hierclust")
result <- melt(result, id = "sd")
ggplot(result, aes(x=sd, y=value, colour=variable)) +
  geom_line(size=2) + 
  geom_point(size=3.5)
ggsave("result.pdf")


# ##################################################
# how many clusters?
# ##################################################
# ConsensusClusterPlus creates a folder which contains the consensus matrices
res_clust<- ConsensusClusterPlus(gene, maxK=6, reps=500, pItem=1, pFeature=1,
                                 clusterAlg="hc",distance="pearson",seed=3,plot="png")



                   
