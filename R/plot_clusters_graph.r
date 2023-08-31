#' Creates a graph plot using the similarity values calculated with ClusterFoldSimilarity().
#'
#' `plotClustersGraph()` Creates a graph plot using the similarity values calculated with ClusterFoldSimilarity().
#'
#' This function will calculate a similarity coeficient using the fold changes of shared genes among clusters of different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster i of sample n and all the target clusters in sample j.
#'
#' @param similarityTable Dataframe. A table obtained from ClusterFoldSimilarity that contains the similarity values as a column "similarity_value" that represents 
#' the similarity of a source cluster to a target cluster.
#' 
#' @return This function plots a graph in which the nodes are clusters from a specific dataset, the edges represent the similarity and the direction
#' of that similarity between clusters.
#' 
#' @examples
#' if (requireNamespace("Seurat") & requireNamespace("SeuratObject")){
#' library(ClusterFoldSimilarity)
#' library(Seurat)
#' library(SeuratObject)
#' # data dimensions
#' nfeatures <- 2000; ncells <- 400
#' # single-cell 1
#' counts <- matrix(runif(nfeatures * ncells, 1, 1e4), nfeatures)
#' rownames(counts) <- paste0("gene",seq(nfeatures))
#' colnames(counts) <- paste0("cell",seq(ncells))
#' colData <- data.frame(cluster=sample(c("Cluster1","Cluster2","Cluster3"),size = ncells,replace = TRUE),
#'                      row.names=paste0("cell",seq(ncells)))
#' seu.1 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
#' Idents(object = seu.1) <- "cluster"
#' # single-cell 2
#' counts <- matrix(runif(nfeatures * ncells, 1, 1e4), nfeatures)
#' rownames(counts) <- paste0("gene",seq(nfeatures))
#' colnames(counts) <- paste0("cell",seq(ncells))
#' colData <- data.frame(cluster=sample(c("Cluster1","Cluster2","Cluster3","Cluster4"),size = ncells,replace = TRUE),
#'                       row.names=paste0("cell",seq(ncells)))
#' seu.2 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
#' Idents(object = seu.2) <- "cluster"
#' # Create a list with the unprocessed single-cell datasets
#' singlecell.object.list <- list(seu.1, seu.2)
#' 
#' similarity.table <- clusterFoldSimilarity(sceList = singlecell.object.list, sampleNames = c("sc1","sc2"))
#' head(similarity.table)
#' plotClustersGraph(similarityTable=similarity.table)
#' }
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom igraph graph_from_data_frame layout_with_fr V
#' @importFrom scales alpha
#' @importFrom graphics legend par
#' @export
plotClustersGraph <- function(similarityTable=NULL){
  ## We will use the similarity matrix to create a graph
  if(is.null(similarityTable) | !is.data.frame(similarityTable)){
    stop("A similarity table from clusterFoldSimilarity was not found.")
  }
  df <- similarityTable
  from <- paste(paste("D", df$datasetL, sep="."), paste0("C", df$clusterL), sep=".")
  to <- paste(paste("D", df$datasetR, sep="."), paste0("C", df$clusterR), sep=".")
  ## igraph
  relations <- data.frame(from=from, to=to, sim=df$similarityValue)
  g <- igraph::graph_from_data_frame(relations, directed=TRUE)
  l <- igraph::layout_with_fr(g)
  ## Retrieve the name and order of the datasets
  setNames <- unique(df$datasetL)
  cl <- as.numeric(apply(table(df$datasetL, df$clusterL)[setNames,] != 0, 1, sum))
  ## Personalized colors  
  clusterColorsPalete <- c("#FDD49E", "#D5BADB", "#7EB6D9", "#DBECDA", "#F28D35", "#4AA147", "#86608E",
                          "#3C7DA6", "#DE77AE", "#D9E8F5", "#92C791", "#D94D1A", "#F2D377")
  clCodes <- clusterColorsPalete[seq(length(setNames))]
  par(mar=c(1, 5, 5, 1));
  base::plot(g, vertex.label=igraph::V(g)$name, edge.arrow.size=.3, layout=l,
             vertex.color=scales::alpha(rep(clCodes, cl), alpha=0.8), edge.color="black",
             vertex.size=14, vertex.frame.color=NA, vertex.label.color="black", 
             vertex.label.cex=0.6, vertex.label.dist=0, edge.curved=0.2)
  # legend('topleft', legend=paste0('Dataset ', unique(df$datasetL)),
  #        fill=clCodes, xpd=TRUE, inset=c(-0.1, -0.1), cex=0.8)
}
