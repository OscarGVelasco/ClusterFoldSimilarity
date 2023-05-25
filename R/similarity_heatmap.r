#' Plot a heatmap of the similarity values obtained using cluster fold similarity
#'
#' `similarityHeatmap()` returns a ggplot heatmap representing the similarity values between pairs of clusters as obtained from \link[ClusterFoldSimilarity]{clusterFoldSimilarity}.
#'
#' This function plots a heatmap using ggplot. It is intended to be used with the output table from \link[ClusterFoldSimilarity]{clusterFoldSimilarity}, which includes the columns: datasetL (the dataset used for comparison)
#' datasetR (the dataset against datasetL has been contrasted), clusterL (clusters from datasetL), clusterR (clusters from datasetR) and the similarityValue.
#'
#' @param similarityTable A DataFrame containing the similarities between all possible pairs of single cell samples obtained with \link[ClusterFoldSimilarity]{clusterFoldSimilarity} using the option n_top=Inf.
#' @param mainDataset Numeric. Specify the main dataset (y axis). It corresponds with the datasetL column from the similarityTable
#' @param otherDatasets Numeric. Specify some specific dataset to be ploted along the mainDataset (x axis, default: all other datasets found on datasetR column from similarity_table).
#' 
#' @return The function returns a heatmap ggplot object. 
#' 
#' @examples
#' if (!requireNamespace("Seurat", quietly=TRUE)){stop("Seurat package not available", call.=FALSE)}
#' if (!requireNamespace("SingleCellExperiment", quietly=TRUE)){stop("SingleCellExperiment package not available", call.=FALSE)}
#' library(ClusterFoldSimilarity)
#' # Mouse brain single-cell RNA-seq 1 from Romanov et. al.
#' mouse.brain.romanov <- scRNAseq::RomanovBrainData(ensembl = TRUE)
#' colnames(mouse.brain.romanov) <- colData(mouse.brain.romanov)$cellID
#' rownames(colData(mouse.brain.romanov)) <- colData(mouse.brain.romanov)$cellID
#' singlecell.1.seurat <- CreateSeuratObject(counts = counts(mouse.brain.romanov),meta.data = as.data.frame(colData(mouse.brain.romanov)))
#' 
#' # Mouse brain single-cell RNA-seq 2 from Zeisel et. al.
#' mouse.brain.zei <- scRNAseq::ZeiselBrainData(ensembl = TRUE)
#' singlecell.2.seurat <- CreateSeuratObject(counts = counts(mouse.brain.zei),meta.data = as.data.frame(colData(mouse.brain.zei)))
#' 
#' # Create a list with the unprocessed single-cell datasets
#' singlecell.object.list <- list(singlecell.1.seurat,singlecell.2.seurat)
#' # Apply the same processing to each dataset and return a list of single-cell analysis
#' singlecell.object.list <- lapply(X = singlecell.object.list, FUN = function(x){
#'   x <- NormalizeData(x)
#'   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
#'   x <- ScaleData(x,features = VariableFeatures(x))
#'   x <- RunPCA(x, features = VariableFeatures(object = x))
#'   x <- FindNeighbors(x, dims = seq(10))
#'   x <- FindClusters(x, resolution = 0.1)
#' })
#' # singlecell.object.list consist of 2 or more single-cell objects
#' # Using topN = Inf by default plots a heatmap using the similarity values: 
#' similarity.table.all <- clusterFoldSimilarity(sceList=singlecell.object.list, topN=Inf)
#' # Using the dataset 2 as a reference on the Y-axis of the heatmap:
#' similarityHeatmap(similarityTable=similarity.table.all, mainDataset=2)
#' 
#' @author Oscar Gonzalez-Velasco
#' @import ggplot2
#' @export
similarityHeatmap <- function(similarityTable = NULL,
                               mainDataset=NULL,
                               otherDatasets=NULL){
myTable <- similarityTable
if(is.null(mainDataset)){
  # If no dataset is specified, we pick the first label from the datasetL values.
  mainDataset <- unique(myTable$datasetL)[1]
}
subMyTable <- myTable[myTable$datasetL %in% c(mainDataset),]
if(!is.null(otherDatasets)){subMyTable <- subMyTable[subMyTable$datasetR %in% c(otherDatasets),]}
# Order the rows by using factor levels:
subMyTable$clusterL <- factor(subMyTable$clusterL, levels=unique(subMyTable$clusterL),ordered = TRUE)
foundDataset <- unique(subMyTable$datasetR)
labels_personaliz <- paste("Dataset",foundDataset)
names(labels_personaliz) <- foundDataset
ylabel <- paste("Dataset",mainDataset,"clusterL")

g <- ggplot2::ggplot(subMyTable, ggplot2::aes(y=clusterL, x=clusterR, fill= similarityValue)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low="#7855ba", # Dark Violet - Low
                                mid="#ffffff", # White - Medium
                                high="#8a0d00", # Dark red - High
                                guide="colorbar",
                                name="Similarity value") + 
  {if(length(foundDataset)>1)ggplot2::facet_grid(~ datasetR, scales="free_x", space="free_x", 
                                                 labeller=ggplot2::as_labeller(labels_personaliz))}+
  {if(length(foundDataset)==1)ggplot2::xlab(paste(labels_personaliz, "clusterR"))}+
  ggplot2::theme_minimal() +
  ggplot2::ylab(ylabel) +
  ggplot2::theme(
    strip.text=ggplot2::element_text(face="bold", size=ggplot2::rel(0.8)),
    strip.background=ggplot2::element_rect(fill="white", colour="black", size=0.6),
    axis.title.x=ggplot2::element_text(),
    axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)
  )
return(g)
}
