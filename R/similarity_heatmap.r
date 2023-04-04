#' Plot a heatmap of the similarity values obtained using cluster fold similarity
#'
#' `similarity_heatmap()` returns a ggplot heatmap representing the similarity values between pairs of clusters as obtained from \link[ClusterFoldSimilarity]{cluster_fold_similarity}.
#'
#' This function plots a heatmap using ggplot. It is intended to be used with the output table from \link[ClusterFoldSimilarity]{cluster_fold_similarity}, which includes the columns: dataset_l (the dataset used for comparison)
#' dataset_r (the dataset against dataset_l has been contrasted), cluster_l (clusters from dataset_l), cluster_r (clusters from dataset_r) and the similarity_value.
#'
#' @param similarity_table A DataFrame containing the similarities between all possible pairs of single cell samples obtained with \link[ClusterFoldSimilarity]{cluster_fold_similarity} using the option n_top=Inf.
#' @param main_dataset Numeric. Specify the main dataset (y axis). It corresponds with the dataset_l column from the similarity_table.
#' @param other_datasets Numeric. Specify some specific dataset to be ploted along the main_dataset (x axis, default: all other datasets found on dataset_r column from similarity_table).
#' 
#' @return The function returns a heatmap ggplot object. 
#' 
#' @examples
#' library(Seurat)
#' library(scRNAseq)
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
#' # Using top_n = Inf by default plots a heatmap using the similarity values: 
#' similarity.table.all <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n = Inf)
#' # Ploting the dataset 2 on the Y-axis:
#' similarity_heatmap(similarity_table = similarity.table.all,main_dataset = 2)
#' 
#' @author Oscar Gonzalez-Velasco
#' @export
similarity_heatmap <- function(similarity_table = NULL,
                               main_dataset=NULL,
                               other_datasets=NULL){
my.table <- similarity_table
if(is.null(main_dataset)){
  # If no dataset is specified, we pick the first label from the dataset_l values.
  main_dataset <- unique(my.table$dataset_l)[1]
}
sub.mt <- my.table[my.table$dataset_l %in% c(main_dataset),]
if(!is.null(other_datasets)){sub.mt <- sub.mt[sub.mt$dataset_r %in% c(other_datasets),]}
# Order the rows by using factor levels:
sub.mt$cluster_l <- factor(sub.mt$cluster_l, levels=unique(sub.mt$cluster_l),ordered = TRUE)
found_dataset <- unique(sub.mt$dataset_r)
labels_personaliz <- paste("Dataset",found_dataset)
names(labels_personaliz) <- found_dataset
ylabel <- paste("Dataset",main_dataset,"cluster_l")

g <- ggplot2::ggplot(sub.mt, ggplot2::aes(y=cluster_l, x=cluster_r, fill= similarity_value)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#7855ba", # Dark Violet - Low
                                mid = "#ffffff", # White - Medium
                                high =   "#8a0d00", # Dark red - High
                                guide="colorbar",
                                name="Similarity value") + 
  {if(length(found_dataset)>1)ggplot2::facet_grid(~ dataset_r, scales = "free_x", space = "free_x",labeller = ggplot2::as_labeller(labels_personaliz))}+
  {if(length(found_dataset)==1)ggplot2::xlab(paste(labels_personaliz,"cluster_r"))}+
  ggplot2::theme_minimal() +
  ggplot2::ylab(ylabel) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.8)),
    strip.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.6),
    axis.title.x=ggplot2::element_text(),
    axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)
  )

return(g)
}
