#' Find cell-group communities by constructing and clustering a directed graph using the similarity values calculated by ClusterFoldSimilarity()
#'
#' `findCommunitiesSimmilarity()` Find communities by constructing and clustering graph using the similarity values calculated by ClusterFoldSimilarity().
#'
#' This function will group together nodes of the network into communities using the InfoMap community detection algorithm.
#'
#' @param similarityTable Dataframe. A table obtained from ClusterFoldSimilarity that contains the similarity values as a column "similarityValue" that represents 
#' the similarity of a source cluster to a target cluster.
#' @param nIterations Numeric. The number of iterations for the Leiden clustering algorithm (default 10).
#' @param leidenResolution Numeric. The resolution parameter for the Leiden clustering algorithm (default 1.2).
#' 
#' @return This function returns a data frame with the community that each node of the network (cell groups defined by the user) belongs to, and plots a graph in 
#' which the nodes are clusters from a specific dataset, the edges represent the similarity and the direction of that similarity between clusters.
#' \tabular{ll}{
#'    \code{sample} \tab The sample name. \cr
#'    \tab \cr
#'    \code{group} \tab The group/cluster from sample defined by the user. \cr
#'    \tab \cr
#'    \code{community} \tab Community group to which the sample:group belongs to.  \cr
#' }
#' @examples
#' if (requireNamespace("Seurat") & requireNamespace("SeuratObject")){
#' library(ClusterFoldSimilarity)
#' library(Seurat)
#' library(SeuratObject)
#' # data dimensions
#' nfeatures <- 2000; ncells <- 400
#' # single-cell 1
#' counts <- matrix(rpois(n=nfeatures * ncells, lambda=10), nfeatures)
#' rownames(counts) <- paste0("gene",seq(nfeatures))
#' colnames(counts) <- paste0("cell",seq(ncells))
#' colData <- data.frame(cluster=sample(c("Cluster1","Cluster2","Cluster3"),size = ncells,replace = TRUE),
#'                      row.names=paste0("cell",seq(ncells)))
#' seu1 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
#' Idents(object = seu1) <- "cluster"
#' # single-cell 2
#' counts <- matrix(rpois(n=nfeatures * ncells, lambda=10), nfeatures)
#' rownames(counts) <- paste0("gene",seq(nfeatures))
#' colnames(counts) <- paste0("cell",seq(ncells))
#' colData <- data.frame(cluster=sample(c("Cluster1","Cluster2","Cluster3","Cluster4"),size = ncells,replace = TRUE),
#'                       row.names=paste0("cell",seq(ncells)))
#' seu2 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
#' Idents(object = seu2) <- "cluster"
#' # Create a list with the unprocessed single-cell datasets
#' singlecellObjectList <- list(seu1, seu2)
#' 
#' similarityTable <- clusterFoldSimilarity(scList = singlecellObjectList, sampleNames = c("sc1","sc2"))
#' head(similarityTable)
#' findCommunitiesSimmilarity(similarityTable=similarityTable)
#' }
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom igraph graph_from_data_frame layout_with_fr V cluster_infomap modularity membership
#' @importFrom scales alpha
#' @importFrom graphics legend par
#' @export
findCommunitiesSimmilarity <- function(similarityTable = NULL, nIterations = 10, leidenResolution = 1) 
{
  if (!is.data.frame(similarityTable)) {
    stop("similarityTable has to be a dataframe return by clusterFoldSimilarity")
  }
  else if (!(c("similarityValue") %in% colnames(similarityTable))) {
    stop("similarityTable has to be a dataframe return by clusterFoldSimilarity")
  }
  df <- similarityTable
  df[which(df$similarityValue < 0), "similarityValue"] <- 1e-05
  from <- paste(df$datasetL, df$clusterL, sep = "_group_")
  to <- paste(df$datasetR, df$clusterR, sep = "_group_")
  relations <- data.frame(from = from, to = to, weight = df$similarityValue)
  g <- igraph::graph_from_data_frame(relations, directed = TRUE)
  gUndirected <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "sum")
  communities <- igraph::cluster_leiden(gUndirected, resolution_parameter = leidenResolution, n_iterations = nIterations)
  comunityMermership <- igraph::membership(communities)
  nCommunities <- length(unique(comunityMermership))
  message("Number of communities: ", nCommunities)
  l <- igraph::layout_with_gem(g)
  setNames <- unique(df$datasetL)
  cl <- as.numeric(apply(table(df$datasetL, df$clusterL)[setNames,] != 0, 1, sum))
  ## Personalized colors  
  clusterColorsPalete <- c("#FDD49E", "#D5BADB", "#7EB6D9", "#DBECDA", "#F28D35", "#4AA147", "#86608E", 
                           "#3C7DA6", "#DE77AE", "#D9E8F5", "#92C791", "#D94D1A", "#F2D377")
  clCodes <- clusterColorsPalete[seq(length(setNames))]
  # Define the colors for the communities:
  groupColores <- hcl.colors(nCommunities, palette = "Set2", alpha = 0.2, rev = FALSE, fixup = TRUE)
  groupColoresBorder <- hcl.colors(nCommunities, palette = "Set2", alpha = 0.6, rev = FALSE, fixup = TRUE)
  par(mar = c(1, 1, 3, 1))
  base::plot(communities, g, layout = l, vertex.label = gsub(pattern = "_group_", replacement = ".", x = igraph::V(g)$name), 
             vertex.color = scales::alpha(rep(clCodes,  cl), alpha = 0.8), edge.color = scales::alpha("black", alpha = 0.6), col = scales::alpha(rep(clCodes, cl), alpha = 0.8), 
             vertex.size = 15, edge.arrow.size = 0.4, vertex.frame.color = NA, 
             vertex.label.color = "black", mark.col = groupColores, 
             mark.border = groupColoresBorder, vertex.label.cex = 0.9, 
             vertex.label.dist = 0, edge.curved = 0.2)
  legend("topleft", horiz = FALSE, y.intersp = 0.8, x.intersp = 0.5, box.lwd = 0,box.col = "white",
         legend = paste0("Dataset ", unique(df$datasetL)), fill = clCodes,
         xpd = TRUE, inset = c(0, 0), cex = 0.8)
  # Return the module info as a data.frame
  samples <- unlist(lapply(strsplit(names(comunityMermership), "_group_"), function(nameEntry) nameEntry[1]))
  groups <- unlist(lapply(strsplit(names(comunityMermership), "_group_"), function(nameEntry) nameEntry[2]))
  dfMemberships <- data.frame(sample = samples, group = groups, community = c(comunityMermership))
  return(dfMemberships)
}
