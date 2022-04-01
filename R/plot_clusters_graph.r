#' Creates a graph plot using the similarity values calculated with ClusterFoldSimilarity().
#'
#' `plot_clusters_graph()` Creates a graph plot using the similarity values calculated with ClusterFoldSimilarity().
#'
#' This function will calculate a similarity coeficient using the fold changes of shared genes among clusters of different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster i of sample n and all the target clusters in sample j.
#'
#' @param similarity.table Dataframe. A table obtained from ClusterFoldSimilarity that contains the similarity values as a column "similarity_value" that represents 
#' the similarity of a source cluster to a target cluster.
#' 
#' @return This function plots a graph in which the nodes are clusters from a specific dataset, the edges represent the similarity and the direction
#' of that similarity between clusters.
#' 
#' @export
plot_clusters_graph <- function(similarity.table = NULL){
  
  # We will use the similarity matrix to create a graph
  if(is.null(similarity.table) | !is.data.frame(similarity.table)){
    stop("A similarity table from clusterFoldSimilarity was not found.")
  }  
  df <- similarity.table
  from <- paste(paste0("D", df$dataset_l), paste0("C", df$cluster_l), sep = ".")
  to <- paste(paste0("D", df$dataset_r), paste0("C", df$cluster_r), sep = ".")
  # igraph
  require(igraph,quietly = T)
  relations <- data.frame(from = from, to = to, sim = df$similarity_value)
  g <- igraph::graph_from_data_frame(relations, directed = T)
  cl <- as.numeric(apply(table(df$dataset_l, df$cluster_l) != 0, 1, sum))
  cluster_colors <- c("#7BCCC4","#FEB24C","#74C476","#FDDBC7","#D0D1E6","#FDBF6F","#D1E5F0","#67A9CF","#E08214","#800026","#006837","#0570B0",
                      "#FC8D62", "#045A8D", "#02818A", "#8C96C6", "#CCEBC5", "#E31A1C", "#78C679", "#A8DDB5", "#FDD49E", "#DE77AE", "#B3B3B3", "#EF6548",
                      "#D73027")
  cl_codes <- sample(cluster_colors,size = length(unique(df$dataset_l)))
  par(mar=c(2,1,1,1));base::plot(g, vertex.label = V(g)$name, edge.arrow.size = .3, 
       vertex.color = rep(cl_codes, cl),edge.color = "black",
       vertex.size = 12, vertex.frame.color = NA, vertex.label.color = "black", 
       vertex.label.cex = 0.6, vertex.label.dist = 0, edge.curved = 0.2);legend('topleft', legend = paste0('Dataset ', 1:length(unique(df$dataset_l))),fill = cl_codes)
}
