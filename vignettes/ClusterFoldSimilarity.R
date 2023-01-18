## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(eval=FALSE,warning=FALSE, error=FALSE, message=FALSE)

## ----construct----------------------------------------------------------------
#  library(Seurat)
#  library(scRNAseq)
#  
#  # Mouse brain single-cell RNA-seq 1 from Romanov et. al.
#  mouse.brain.romanov <- scRNAseq::RomanovBrainData(ensembl = T,location = F)
#  colnames(mouse.brain.romanov) <- colData(mouse.brain.romanov)$cellID
#  rownames(colData(mouse.brain.romanov)) <- colData(mouse.brain.romanov)$cellID
#  singlecell.1.seurat <- CreateSeuratObject(counts = counts(mouse.brain.romanov),meta.data = as.data.frame(colData(mouse.brain.romanov)))
#  
#  # Mouse brain single-cell RNA-seq 2 from Zeisel et. al.
#  mouse.brain.zei <- scRNAseq::ZeiselBrainData(ensembl = T,location = F)
#  singlecell.2.seurat <- CreateSeuratObject(counts = counts(mouse.brain.zei),meta.data = as.data.frame(colData(mouse.brain.zei)))

## -----------------------------------------------------------------------------
#  # Create a list with the unprocessed single-cell datasets
#  singlecell.object.list <- list(singlecell.1.seurat,singlecell.2.seurat)
#  # Apply the same processing to each dataset and return a list of single-cell analysis
#  singlecell.object.list <- lapply(X = singlecell.object.list, FUN = function(x){
#    x <- NormalizeData(x)
#    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
#    x <- ScaleData(x,features = VariableFeatures(x))
#    x <- RunPCA(x, features = VariableFeatures(object = x))
#    x <- FindNeighbors(x, dims = 1:10)
#    x <- FindClusters(x, resolution = 0.1)
#  })

## -----------------------------------------------------------------------------
#  library(ClusterFoldSimilarity)
#  similarity.table <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n = 1)
#  head(similarity.table)

## -----------------------------------------------------------------------------
#  label1 = "level1.class" # name of label of data 1
#  label2 = "level1class" # name of label of data 2
#  
#  apply(similarity.table[similarity.table$dataset_l == 1,],1,function(x){
#    n1 = names(which.max(table(singlecell.object.list[[as.numeric(x["dataset_l"])]]@meta.data[singlecell.object.list[[as.numeric(x["dataset_l"])]]@meta.data$seurat_clusters == x["cluster_l"],label1])))
#    n2 = names(which.max(table(singlecell.object.list[[as.numeric(x["dataset_r"])]]@meta.data[singlecell.object.list[[as.numeric(x["dataset_r"])]]@meta.data$seurat_clusters == x["cluster_r"],label2])))
#    return(paste("dataset 1 cluster",x["cluster_l"],"top cell.type:",n1,"VS dataset 2 cluster",x["cluster_r"],"top cell.type:",n2))
#    })

## -----------------------------------------------------------------------------
#  # Retrieve the top 3 similar cluster for each of the clusters:
#  similarity.table.3top <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n = 3)
#  head(similarity.table.3top)

## -----------------------------------------------------------------------------
#  # Retrieve the top 5 features that contribute the most to the similarity between each pair of clusters:
#  similarity.table.5top.features <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n_genes = 5)
#  head(similarity.table.5top.features, n=10)

## -----------------------------------------------------------------------------
#  similarity.table.all.values <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n = Inf)
#  dim(similarity.table.all.values)

## -----------------------------------------------------------------------------
#  dataset1= 1
#  dataset2= 2
#  similarity.table.2 <- similarity.table.all.values %>%
#                        filter(dataset_l == dataset1 & dataset_r == dataset2) %>%
#                        arrange(desc(as.numeric(cluster_l)),as.numeric(cluster_r))
#  cls <- unique(similarity.table.2$cluster_l)
#  cls2 <- unique(similarity.table.2$cluster_r)
#  similarity.matrix.all <- t(matrix(similarity.table.2$similarity_value,ncol=length(unique(similarity.table.2$cluster_l))))
#  rownames(similarity.matrix.all) <- cls
#  colnames(similarity.matrix.all) <- cls2
#  similarity.matrix.all

## ----setup--------------------------------------------------------------------
#  library(ClusterFoldSimilarity)

