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
#' @param highlightTop Boolean. If the top 2 similarity values should be highlighted on the heatmap (default: TRUE)
#' 
#' @return The function returns a heatmap ggplot object. 
#' 
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
#' # Using topN = Inf by default plots a heatmap using the similarity values: 
#' similarityTableAll <- clusterFoldSimilarity(scList=singlecellObjectList, topN=Inf)
#' # Using the dataset 2 as a reference on the Y-axis of the heatmap:
#' similarityHeatmap(similarityTable=similarityTableAll, mainDataset=2, highlightTop=FALSE)
#' }
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom dplyr %>% filter arrange
#' @importFrom ggdendro ggdendrogram
#' @import ggplot2
#' @importFrom cowplot align_plots plot_grid
#' @importFrom reshape2 melt
#' @importFrom stats hclust dist setNames
#' @export
similarityHeatmap <- function(similarityTable=NULL,
                              mainDataset=NULL,
                              otherDatasets=NULL,
                              highlightTop=TRUE){
  myTable <- similarityTable
  if(is.null(mainDataset)){
    # If no dataset is specified, we pick the first label from the datasetL values.
    mainDataset <- unique(myTable$datasetL)[1]
  }
  subMyTable <- myTable[myTable$datasetL %in% c(mainDataset),]
  if(!is.null(otherDatasets)){subMyTable <- subMyTable[subMyTable$datasetR %in% c(otherDatasets),]}
  # Order the rows by using factor levels:
  subMyTable$clusterL <- factor(subMyTable$clusterL, levels=unique(subMyTable$clusterL), ordered=TRUE)
  foundDataset <- unique(subMyTable$datasetR)
  labels_personaliz <- paste("Dataset", foundDataset)
  names(labels_personaliz) <- foundDataset
  ylabel <- paste("Dataset", mainDataset, "clusterL")
  # Create the dendrogram for Rows:
  nDistinctClasses <- length(unique(paste0(subMyTable[,"datasetR"], subMyTable[,"clusterR"])))
  subMyTable <- subMyTable[order(subMyTable$clusterR),]
  rows <- unique(subMyTable[,"clusterL"])
  similarityMatrixAll <- matrix(subMyTable[,"similarityValue"], ncol=nDistinctClasses)
  rownames(similarityMatrixAll) <- rows
  dendoGlobalRows <- stats::hclust(stats::dist(similarityMatrixAll))
  orderingHierarchicalRows <- dendoGlobalRows$labels[dendoGlobalRows$order]
  # dendrogramRows <- ggdendro::ggdendrogram(dendoGlobalRows, labels=FALSE, rotate=TRUE) +  
  #   ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 0), units="npc"),
  #         axis.text.y=ggplot2::element_blank(),
  #         axis.text.x=ggplot2::element_blank())
  # We split the similarity table by reference datasets:
  similarityDfList <- split(subMyTable, f=subMyTable$datasetR )
  # Create one plot per reference dataset:
  listOfPlots <- lapply(similarityDfList, function(df){
    # Contruct the Dendogram:
    dataset2 <- unique(df$datasetR)
    dataset1 <- mainDataset
    similarityMatrix <- df %>% 
      dplyr::filter(datasetL == dataset1 & datasetR == dataset2) %>% 
      dplyr::arrange(desc(clusterL), clusterR)
    
    similarityMatrix <- similarityMatrix[order(similarityMatrix$clusterR),]
    
    rows <- unique(similarityMatrix$clusterL)
    columns <- as.character(unique(similarityMatrix$clusterR))
    similarityMatrixAll <- matrix(similarityMatrix$similarityValue, ncol=length(columns))
    rownames(similarityMatrixAll) <- rows
    colnames(similarityMatrixAll) <- columns
    # Build the dendrogram
    dendo <- stats::hclust(stats::dist(t(similarityMatrixAll)))
    orderingHierarchical <- dendo$labels[dendo$order]
    dendrogram <- ggdendro::ggdendrogram(dendo, labels=FALSE) +  
      ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 0), units="npc"),
            axis.text.y=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank())
    # Matrix to plot:
    df <- reshape2::melt(similarityMatrixAll)
    colnames(df) <- c("clusterL", "clusterR", "similarityValue")
    df <- df %>% dplyr::arrange(factor(clusterR, levels = orderingHierarchical))
    df$clusterR <- factor(df$clusterR, levels = orderingHierarchical)
    df$clusterL <- factor(df$clusterL, levels = orderingHierarchicalRows)
    if(isTRUE(highlightTop)){
      highlightTop <- data.frame(clusterR=max.col(similarityMatrixAll[orderingHierarchicalRows,orderingHierarchical]),
                                 clusterL=1:nrow(similarityMatrixAll))
      highlightSecond <- data.frame(clusterR=apply(similarityMatrixAll[orderingHierarchicalRows,orderingHierarchical],1,function(row){order(row, decreasing=TRUE)[2]}),
                                    clusterL=1:nrow(similarityMatrixAll))
    }
    if(isFALSE(highlightTop)){
      highlightTop <- stats::setNames(data.frame(matrix(ncol=2, nrow=0)), c("clusterR", "clusterL"))
      highlightSecond <- stats::setNames(data.frame(matrix(ncol=2, nrow=0)), c("clusterR", "clusterL"))
    }
    # Build the plot
    ggHeatmap <- ggplot2::ggplot(df) +
      ggplot2::aes(y=clusterL, x=clusterR, fill=similarityValue) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low="#7855ba", # Dark Violet - Low
                                    mid="#ffffff", # White - Medium
                                    high="#8a0d00", # Dark red - High
                                    guide="colorbar",
                                    name="Similarity value") +
      ggplot2::geom_rect(data=highlightSecond, size=0.6, fill=NA, colour="#990B00",
                         ggplot2::aes(xmin=clusterR - 0.5, xmax=clusterR + 0.5, ymin=clusterL - 0.5, ymax=clusterL + 0.5)) +
      ggplot2::geom_rect(data=highlightTop, size=0.6, fill=NA, colour="black",
                         ggplot2::aes(xmin=clusterR - 0.5, xmax=clusterR + 0.5, ymin=clusterL - 0.5, ymax=clusterL + 0.5)) +
      ggplot2::theme_minimal() +
      ggplot2::ylab(ylabel) +
      ggplot2::xlab(paste("Dataset",dataset2,"clusterR")) +
      ggplot2::theme(
        strip.text=ggplot2::element_text(face="bold", size=ggplot2::rel(0.8)),
        strip.background=ggplot2::element_rect(fill="white", colour="black", linewidth=0.6),
        axis.title.x=ggplot2::element_text(),
        axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1),
        legend.position="none", 
        plot.margin=ggplot2::unit(c(0, 0, 0, 0), units="npc")
      )
    return(list(dendrogram, ggHeatmap))
  })
  if(length(foundDataset) > 1){
    heatmaps <- cowplot::align_plots(plotlist=lapply(listOfPlots, function(plot){plot[[2]]}), align="h", axis="b")
    dendros <- cowplot::align_plots(plotlist=lapply(listOfPlots, function(plot){plot[[1]]}), align="h", axis="b")
    finalPlot <- cowplot::plot_grid(plotlist=c(dendros, heatmaps), align="v", axis="tb", rel_heights=c(0.2, 1))
  } else{
    finalPlot <- cowplot::plot_grid(listOfPlots[[1]][[1]], listOfPlots[[1]][[2]], nrow = 2, align="v", axis="tb", rel_heights=c(0.2, 1))
  }
    return(finalPlot)
  }
