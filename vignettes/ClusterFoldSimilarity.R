## ----options, include=FALSE, echo=FALSE---------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(eval=TRUE,warning=FALSE, error=FALSE, message=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("OscarGVelasco/ClusterFoldSimilarity")

## ----setup--------------------------------------------------------------------
library(ClusterFoldSimilarity)

## ----construct----------------------------------------------------------------
library(Seurat)
library(scRNAseq)
library(dplyr)
# Human pancreatic single cell data 1
pancreas_muraro <- scRNAseq::MuraroPancreasData(ensembl = FALSE)
table(colData(pancreas_muraro)$label); dim(pancreas_muraro)
pancreas_muraro <- pancreas_muraro[,rownames(colData(pancreas_muraro)[!is.na(colData(pancreas_muraro)$label),])]
colData(pancreas_muraro)$cell.type <- colData(pancreas_muraro)$label
rownames(pancreas_muraro) <- unlist(lapply(strsplit(rownames(pancreas_muraro),split = "__"),function(x)x[[1]]))
singlecell.1.seurat <- CreateSeuratObject(counts = counts(pancreas_muraro),meta.data = as.data.frame(colData(pancreas_muraro)))

# Human pancreatic single cell data 2
pancreas_baron <- scRNAseq::BaronPancreasData(which = "human",ensembl = FALSE)
table(colData(pancreas_baron)$label); dim(pancreas_baron)
colData(pancreas_baron)$cell.type <- colData(pancreas_baron)$label
singlecell.2.seurat <- CreateSeuratObject(counts = counts(pancreas_baron),meta.data = as.data.frame(colData(pancreas_baron)))


## -----------------------------------------------------------------------------
# Create a list with the unprocessed single-cell datasets
singlecell.object.list <- list(singlecell.1.seurat,singlecell.2.seurat)
# Apply the same processing to each dataset and return a list of single-cell analysis
singlecell.object.list <- lapply(X = singlecell.object.list, FUN = function(x){
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x <- ScaleData(x,features = VariableFeatures(x))
x <- RunPCA(x, features = VariableFeatures(object = x))
x <- FindNeighbors(x, dims = seq(16))
x <- FindClusters(x, resolution = 0.4)
})

## -----------------------------------------------------------------------------
# Assign cell-type annotated from the original study to the cell labels:
Idents(singlecell.object.list[[1]]) <- factor(singlecell.object.list[[1]][[]][,"cell.type"])

library(ClusterFoldSimilarity)
similarity.table <- clusterFoldSimilarity(sceList=singlecell.object.list, sampleNames = c("human","human.NA"),
                                          topN=1, nSubsampling = 24)
head(similarity.table)

## -----------------------------------------------------------------------------

type.count <- singlecell.object.list[[2]][[]] %>% 
  group_by(seurat_clusters) %>% 
  count(cell.type) %>%
  arrange(desc(n), .by_group = TRUE) %>% 
  filter(row_number()==1)

cbind.data.frame(type.count, 
                 matched.type = rep(table(type.count$seurat_clusters), x = similarity.table[similarity.table$datasetL == "human.NA",]$clusterR))


## -----------------------------------------------------------------------------
# Retrieve the top 3 similar cluster for each of the clusters:
similarity.table.3top <- clusterFoldSimilarity(sceList=singlecell.object.list, topN=3,
                                             sampleNames = c("human","human.NA"), nSubsampling = 24)
head(similarity.table.3top)

## -----------------------------------------------------------------------------
# Retrieve the top 5 features that contribute the most to the similarity between each pair of clusters:
similarity.table.5top.features <- clusterFoldSimilarity(sceList=singlecell.object.list, topNFeatures=5, nSubsampling = 24)
head(similarity.table.5top.features, n=10)

## -----------------------------------------------------------------------------
similarity.table.all.values <- clusterFoldSimilarity(sceList=singlecell.object.list, sampleNames = c("human","human.NA"), topN=Inf)
dim(similarity.table.all.values)

## -----------------------------------------------------------------------------
library(dplyr)
dataset1 <- "human"
dataset2 <- "human.NA"
similarity.table.2 <- similarity.table.all.values %>% 
filter(datasetL == dataset1 & datasetR == dataset2) %>% 
arrange(desc(as.numeric(clusterL)), as.numeric(clusterR))
cls <- unique(similarity.table.2$clusterL)
cls2 <- unique(similarity.table.2$clusterR)
similarity.matrix.all <- t(matrix(similarity.table.2$similarityValue, ncol=length(unique(similarity.table.2$clusterL))))
rownames(similarity.matrix.all) <- cls
colnames(similarity.matrix.all) <- cls2
similarity.matrix.all

## -----------------------------------------------------------------------------
# Mouse pancreatic single cell data
pancreas_baron_mm <- scRNAseq::BaronPancreasData(which = "mouse",ensembl = FALSE)
table(colData(pancreas_baron_mm)$label); dim(pancreas_baron_mm)
colData(pancreas_baron_mm)$cell.type <- colData(pancreas_baron_mm)$label
# Translate mouse gene ids to human ids
# *for the sake of simplicity we are going to transform to uppercase all mouse gene names
rownames(pancreas_baron_mm) <- toupper(rownames(pancreas_baron_mm))
# Create seurat object
singlecell.3.seurat <- CreateSeuratObject(counts = counts(pancreas_baron_mm),meta.data = as.data.frame(colData(pancreas_baron_mm)))

# We append the single-cell object to our list
singlecell.object.list[[3]] <- singlecell.3.seurat


## -----------------------------------------------------------------------------

x <- singlecell.object.list[[3]]
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x <- ScaleData(x,features = VariableFeatures(x))
x <- RunPCA(x, features = VariableFeatures(object = x))
x <- FindNeighbors(x, dims = seq(16))
x <- FindClusters(x, resolution = 0.4)
singlecell.object.list[[3]] <- x

# We use the cell labels as a second reference, but we can also use the cluster labels if our interest is to match clusters
Idents(singlecell.object.list[[3]]) <- factor(singlecell.object.list[[3]][[]][,"cell.type"])

# We subset the most variable genes in each experiment
singlecell.object.list.variable <- lapply(singlecell.object.list, function(x){x[VariableFeatures(x),]})

similarity.table.human.mouse <- clusterFoldSimilarity(sceList = singlecell.object.list.variable,
                                                        sampleNames = c("human","human.NA","mouse"),
                                                        topN = 1, 
                                                        nSubsampling = 24,
                                                        parallel = TRUE)

## -----------------------------------------------------------------------------
similarity.table.human.mouse.all <- clusterFoldSimilarity(sceList = singlecell.object.list.variable,
                                                          sampleNames = c("human","human.NA","mouse"),
                                                          topN = Inf, 
                                                          nSubsampling = 24,
                                                          parallel = TRUE)
# We can select which dataset to plot in the Y-axis:
ClusterFoldSimilarity::similarityHeatmap(similarityTable=similarity.table.human.mouse.all, mainDataset="human.NA")

## -----------------------------------------------------------------------------
sessionInfo()

