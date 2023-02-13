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
#' @examples 
#' # singlecell.object.list consist of 2 or more single-cell objects
#' # By default plots a graph using the similarity values: 
#' similarity.table.top <- cluster_fold_similarity(sce_list = singlecell.object.list,top_n = 1)
#' # The same plot can be created using:
#' plot_clusters_graph(similarity.table = similarity.table.top)
#' 
#' @author Oscar Gonzalez-Velasco
#' @export
plot_clusters_graph <- function(similarity.table = NULL){

  # require(igraph)
  # We will use the similarity matrix to create a graph
  if(is.null(similarity.table) | !is.data.frame(similarity.table)){
    stop("A similarity table from clusterFoldSimilarity was not found.")
  }
  df <- similarity.table
  from <- paste(paste("D", df$dataset_l,sep = "."), paste0("C", df$cluster_l), sep = ".")
  to <- paste(paste("D", df$dataset_r,sep = "."), paste0("C", df$cluster_r), sep = ".")
  # igraph
  relations <- data.frame(from = from, to = to, sim = df$similarity_value)
  g <- igraph::graph_from_data_frame(relations, directed = TRUE)
  cl <- as.numeric(apply(table(df$dataset_l, df$cluster_l) != 0, 1, sum))
  cluster_colors <- c("#FDD49E", "#D5BADB","#7EB6D9","#DBECDA","#F28D35","#4AA147","#86608E","#3C7DA6","#DE77AE","#D9E8F5","#92C791",
              "#D94D1A","#F2D377")
  cl_codes <- cluster_colors[seq(length(unique(df$dataset_l)))]
  par(mar=c(1,5,5,1));
  base::plot(g, vertex.label = igraph::V(g)$name, edge.arrow.size = .3, 
       vertex.color = rep(cl_codes, cl),edge.color = "black",
       vertex.size = 14, vertex.frame.color = NA, vertex.label.color = "black", 
       vertex.label.cex = 0.6, vertex.label.dist = 0, edge.curved = 0.2);
  legend('topleft', legend = paste0('Dataset ', unique(df$dataset_l)),
         fill = cl_codes, xpd=TRUE, inset=c(-0.1,-0.1),cex=0.8)
}