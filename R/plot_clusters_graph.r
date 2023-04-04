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
  # Retrieve the name and order of the datasets
  set_names <- unique(df$dataset_l)
  cl <- as.numeric(apply(table(df$dataset_l, df$cluster_l)[set_names,] != 0, 1, sum))
  # Personalized colors  
  cluster_colors <- c("#FDD49E", "#D5BADB","#7EB6D9","#DBECDA","#F28D35","#4AA147","#86608E","#3C7DA6","#DE77AE","#D9E8F5","#92C791",
              "#D94D1A","#F2D377")
  cl_codes <- cluster_colors[seq(length(set_names))]
  par(mar=c(1,5,5,1));
  base::plot(g, vertex.label = igraph::V(g)$name, edge.arrow.size = .3, 
             vertex.color = scales::alpha(rep(cl_codes, cl),alpha = 0.8), edge.color = "black",
             vertex.size = 14, vertex.frame.color = NA, vertex.label.color = "black", 
             vertex.label.cex = 0.6, vertex.label.dist = 0, edge.curved = 0.2);
  legend('topleft', legend = paste0('Dataset ', unique(df$dataset_l)),
         fill = cl_codes, xpd=TRUE, inset=c(-0.1,-0.1),cex=0.8)
}