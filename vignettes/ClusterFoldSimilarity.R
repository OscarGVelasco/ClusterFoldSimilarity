## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(eval=FALSE,warning=FALSE, error=FALSE, message=FALSE)

## ----construct----------------------------------------------------------------
#  library(scran)
#  set.seed(100)
#  sc.experiment.1 <- mouse.brain.zei
#  lib.sce <- librarySizeFactors(sc.experiment.1)
#  # We choose the deconvolutional factor normalization method, we use the computed clusters to calculate the sizefactors
#  clust.sce <- quickCluster(sc.experiment.1) # Normalization by deconvolution to account for different cell type variability
#  sc.experiment.1 <- computeSumFactors(sc.experiment.1, cluster=clust.sce, min.mean=0.1)
#  sc.experiment.1 <- logNormCounts(sc.experiment.1) # Calculate the log10 normalized matrix
#  # We obtain the top 1000 variant genes
#  hvg.top <- getTopHVGs(sc.experiment.1, n=1000)
#  # We obtain the PCA for downstream analysis:
#  sc.experiment.1 <- runPCA(sc.experiment.1)
#  # Clustering
#  library(bluster)
#  kgraph.clusters <- clusterRows(reducedDim(sc.experiment.1, "PCA"),
#                                 TwoStepParam(
#                                   first=KmeansParam(centers=100),
#                                   second=NNGraphParam(k=4)
#                                 )
#  )
#  table(kgraph.clusters)
#  # We set the cluster information for each cell:
#  colLabels(sc.experiment.1) <- factor(kgraph.clusters)
#  ##
#  
#  my.sc.experiment.list = list()
#  my.sc.experiment.list[[1]] <- sc.experiment.1
#  my.sc.experiment.list[[2]] <- sc.experiment.2

## -----------------------------------------------------------------------------
#  library(ClusterFoldSimilarity)
#  similarity.table <- cluster_fold_similarity(sce_list = my.sc.experiment.list,top_n = 1)
#  head(similarity.table)

## -----------------------------------------------------------------------------
#  library(ClusterFoldSimilarity)
#  similarity.table <- cluster_fold_similarity(sce_list = my.sc.experiment.list,top_n = Inf)
#  dim(similarity.table)

## -----------------------------------------------------------------------------
#  dataset1= 1
#  dataset2= 2
#  similarity.table.2 <- similarity.table %>% filter(dataset_l == dataset1 & dataset_r == dataset2) %>% arrange(desc(as.numeric(cluster_l)),as.numeric(cluster_r))
#  cls <- unique(similarity.table.2$cluster_l)
#  cls2 <- unique(similarity.table.2$cluster_r)
#  similarity.table.all <- t(matrix(similarity.table.2$similarity_value,ncol=length(unique(similarity.table.2$cluster_l))))
#  rownames(similarity.table.all) <- cls
#  colnames(similarity.table.all) <- cls2
#  similarity.table.all

## ----setup--------------------------------------------------------------------
#  library(ClusterFoldSimilarity)

