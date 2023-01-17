#' Plot a heatmap of the similarity values obtained using cluster fold similarity
#'
#' `similarity_heatmap()` returns a ggplot heatmap representing the similarity values between pairs of clusters as obtained from \link[ClusterFoldSimilarity]{cluster_fold_similarity}.
#'
#' This function plots a heatmap using ggplot. It is intended to be used with the output table from \link[ClusterFoldSimilarity]{cluster_fold_similarity}, which includes the columns: dataset_l (the dataset used for comparison)
#' dataset_r (the dataset against dataset_l has been contrasted), cluster_l (clusters from dataset_l), cluster_r (clusters from dataset_r) and the similarity_value.
#'
#' @param similarity_table \linkS4class{DataFrame} containing the similarities between all possible pairs of single cell samples obtained with \link[ClusterFoldSimilarity]{cluster_fold_similarity} using the option n_top=Inf.
#' @param main_dataset Numeric. Specify the main dataset (y axis). It corresponds with the dataset_l column from the similarity_table.
#' @param other_datasets Numeric. Specify some specific dataset to be ploted along the main_dataset (x axis, default: all other datasets found on dataset_r column from similarity_table).
#' 
#' @return The function returns a heatmap ggplot object. 
#' 
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
                                guide="colorbar") + 
  {if(length(found_dataset)>1)ggplot2::facet_grid(~ dataset_r, scales = "free_x", space = "free_x",labeller = ggplot2::as_labeller(labels_personaliz))}+
  ggplot2::theme_minimal() +
  ggplot2::ylab(ylabel) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.8)),
    strip.background = ggplot2::element_rect(fill = "white", colour = "black", size = 0.6),
    axis.title.x=ggplot2::element_text()
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.title.y=element_blank()
  )
return(g)
}
