#' Calculate cluster similarity between clusters from different single cell samples.
#'
#' `clusterFoldSimilarity()` returns a dataframe containing the best top similarities between all possible pairs of single cell samples.
#'
#' This function will calculate a similarity coeficient using the fold changes of shared features (e.g.:genes for a single-cell RNA-Seq, peaks for ATAC-Seq) among clusters, or user-defined-groups, from different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster/group i from sample n and all the target clusters/groups in sample j.
#'
#' @param scList List. A list of Single Cell Experiments or Seurat objects. At least 2 are needed. The objects are expected to have cluster or label groups set as identity class.
#' @param sampleNames Character Vector. Specify the sample names, if not a number corresponding with its position on (scList).
#' @param topN Numeric. Specifies the number of target clusters with best similarity to report for each cluster comparison (default 1). If set to Inf, then all similarity values from all possible pairs of clusters are returned.
#' @param topNFeatures Numeric. Number of top features that explains the clusters similarity to report for each cluster comparison (default 1). If topN = Inf then topNFeatures is automatically set to 1.
#' @param nSubsampling Numeric. Number of random sampling of cells to achieve fold change stability (default 15).
#' @param parallel Boolean. Whether to use parallel computing using BiocParallel or not (default FALSE).
#' 
#' @return The function returns a DataFrame containing the best top similarities between all possible pairs of single cell samples. Column values are:
#' \tabular{ll}{
#'    \code{similarityValue} \tab The top similarity value calculated between datasetL:clusterL and datasetR. \cr
#'    \tab \cr
#'    \code{sem} \tab Standar Error of the Mean (SEM) of the scalar contribution coefficients computed for all features. \cr
#'    \tab \cr
#'    \code{w} \tab Weight associated with the similarity score value. \cr
#'    \tab \cr
#'    \code{datasetL} \tab Dataset left, the dataset/sample which has been used to be compared.  \cr
#'    \tab \cr
#'    \code{clusterL} \tab Cluster left, the cluster source from datasetL which has been compared. \cr
#'    \tab \cr
#'    \code{datasetR} \tab Dataset right, the dataset/sample used for comparison against datasetL. \cr
#'    \tab \cr
#'    \code{clusterR} \tab Cluster right, the cluster target from datasetR which is being compared with the clusterL from datasetL. \cr
#'    \tab \cr
#'    \code{topFeatureConserved} \tab The features (e.g.: genes, peaks...) that most contributed to the similarity between clusterL & clusterR. \cr
#'    \tab \cr
#'    \code{featureScore} \tab The similarity score contribution for the specific topFeatureConserved (e.g.: genes, peaks...). \cr
#' }
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
#' counts <- matrix(rpois(n=nfeatures * ncells, lambda=20), nfeatures)
#' rownames(counts) <- paste0("gene",seq(nfeatures))
#' colnames(counts) <- paste0("cell",seq(ncells))
#' colData <- data.frame(cluster=sample(c("Cluster1","Cluster2","Cluster3","Cluster4"),size = ncells,replace = TRUE),
#'                       row.names=paste0("cell",seq(ncells)))
#' seu2 <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = colData)
#' Idents(object = seu2) <- "cluster"
#' # Create a list with the unprocessed single-cell datasets
#' singlecellObjectList <- list(seu1, seu2)
#' 
#' similarityTable <- clusterFoldSimilarity(scList=singlecellObjectList, sampleNames = c("sc1","sc2"))
#' head(similarityTable)
#' }
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix Matrix
#' @export
clusterFoldSimilarity <- function(scList=NULL, sampleNames=NULL, topN=1, topNFeatures=1, nSubsampling=15, parallel=FALSE){
  if(is.null(scList) | (length(scList)<2)){
    stop("At least two Single Cell Experiments are needed for cluster comparison.")
  }
  if(is.infinite(topN)){
    ## if topN clusters to report is set to Inf, then we return ALL similarity values from ALL the possible combination of clusters
    message("Returning similarity values from ALL pairs of clusters")
    message("  *(Note: by using topN=Inf the topNFeatures is set automatically to 1)")
    topNFeatures <- 1
  }
  isSeurat <- FALSE
  isSce <- FALSE
  ## Function starts by loading dependencies
  if(is(scList[[1]], 'Seurat')){
    if (!requireNamespace("Seurat", quietly=TRUE)) {
      stop("Package \"Seurat\" needed for this function to work. Please install it.",
           call.=FALSE)}
    isSeurat <- TRUE
    isSce <- FALSE
  }
  if(is(scList[[1]], 'SingleCellExperiment')){
    if (!requireNamespace("SingleCellExperiment", quietly=TRUE)) {
      stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.",
           call.=FALSE)}
    isSeurat <- FALSE
    isSce <- TRUE
  }
  if (!isSeurat & !isSce){
    stop("One or more objects in the input list is neither of class Seurat nor SingleDataExperiment.")
  }
  functToApply <- base::lapply
  if(isTRUE(parallel)){
    functToApply <- BiocParallel::bplapply
  }
  ## Final var names
  summaryResults <- data.frame(similarityValue=integer(),
                                sem=integer(),
                                w=integer(),
                                datasetL=integer(), 
                                clusterL=integer(),
                                datasetR=integer(),
                                clusterR=integer(),
                                topFeatureConserved=character(),
                                featureScore=integer(),
                                stringsAsFactors=FALSE)
  ## Select common genes in all samples 
  features <- Reduce(intersect,lapply(scList,function(x){rownames(x)}))
  if(length(features) == 0){
    stop("No common features between datasets. Please, select a subset of common variable features among all datasets.")
  }
  if(length(features) < 10){
    warning("The number of common features among datasets is: ",length(features),". More than 50 common features among all datasets is recomended.")
  }
  ## Control filtered level factors using droplevels
  for (nScObject in length(scList)){
    if(isSce){
      if(is.null(SingleCellExperiment::colLabels(scList[[nScObject]]))){
        stop("No clusters specified on colLabels on sample ",nScObject)
      }
      if(is.factor(SingleCellExperiment::colLabels(scList[[nScObject]]))){
        SingleCellExperiment::colLabels(scList[[nScObject]]) <- droplevels(SingleCellExperiment::colLabels(scList[[nScObject]]))
      }else{
        SingleCellExperiment::colLabels(scList[[nScObject]]) <- factor(SingleCellExperiment::colLabels(scList[[nScObject]]))
      }
      }
    if(isSeurat){
      if(is.factor(Seurat::Idents(scList[[nScObject]]))){
        Seurat::Idents(scList[[nScObject]]) <- droplevels(Seurat::Idents(scList[[nScObject]]))
      }else{
        Seurat::Idents(scList[[nScObject]]) <- factor(Seurat::Idents(scList[[nScObject]]))
      }
      }
  }
  spaces <- "                  " ## Empty spaces to clear the append of message report
  ## Save the cluster name identification given by the user
  if(isSce){clusterNames <- lapply(scList,function(x)levels(SingleCellExperiment::colLabels(x)))}
  if(isSeurat){clusterNames <- lapply(scList,function(x)levels(Seurat::Idents(x)))}
  ## Check for dataset names:
  if(is.null(sampleNames)){
    sampleNames <- seq(length(scList))
  }else if(!(length(sampleNames) == length(scList))){
    stop("Number of names given in sampleNames does not match the length of the single-cell experiment list.")
  }
  ## Calculate cluster FoldChange pairwise values:
  markerScList <- list()
  message("Using a common set of", length(features), "features.")
  # Testing the number of draws/subsamplings of cells needed to see *all cells*
  cellDraw <- function(n){
    percentage <- (1/3) * n
    return(sum(rep(1, n+1) / seq(1, n+1)) * (n/percentage))
  }
  if(isSce){
    nOfDraws <- unlist(lapply(seq_len(length(scList)),function(i)lapply(c(table(SingleCellExperiment::colLabels(scList[[i]]))),cellDraw)))
  }
  if(isSeurat){
    nOfDraws <- unlist(lapply(seq_len(length(scList)),function(i)lapply(c(table(Seurat::Idents(scList[[i]]))),cellDraw)))
  }
  # We select as optimal the subsampling n percentile 85% across all cell groups
  message("Using a cell subsampling of n=",nSubsampling," (recomended n=",round(quantile(nOfDraws, probs=0.85)),")")
  message("Computing fold changes.")
  if(isSce){
    markerScList <- functToApply(scList, function(x){
      pairwiseClusterFoldChange(countData=Matrix::Matrix(data=SingleCellExperiment::counts(x), sparse=TRUE), clusters=SingleCellExperiment::colLabels(x), nSubsampling=nSubsampling, functToApply=functToApply) # Raw counts
    })
  }
  if(isSeurat){
    markerScList <- functToApply(scList, function(x){
      pairwiseClusterFoldChange(countData=Matrix::Matrix(data=Seurat::GetAssayData(x, slot="counts"), sparse=TRUE), clusters=Seurat::Idents(x), nSubsampling=nSubsampling, functToApply=functToApply) # Raw counts from Seurat object
    })
  }
  ## Main loop - samples
  intermediateResults <- lapply(seq_len(length(markerScList)),function(sampleN){
    ## Pick a sample in ascending order
    rootSampleList <- markerScList[[sampleN]]
    intermediateResults <- functToApply(seq_along(rootSampleList),function(clusterN){ # Parallel computing if selected
      ## We choose a root cluster and compare it with clusters of samples sampleN+1+...+n
      root <- rootSampleList[[clusterN]]
      for (sampleNTarget in seq(from=1, to=length(markerScList), by=1)){
        if(sampleNTarget == sampleN)next() # If root and target samples are the same, we go for the next sample and skip this loop
        sceComparative <- markerScList[[sampleNTarget]]
        results <- data.frame(similarityValue=integer(),
                              sem=integer(),
                              w=integer(),
                              datasetL=integer(), 
                              clusterL=integer(),
                              datasetR=integer(),
                              clusterR=integer(),
                              topFeatureConserved=character(),
                              featureScore=integer(),
                              stringsAsFactors=FALSE)
          for(clusterNTarget in seq_along(sceComparative)){
            ## The sample for comparing will be the sample_i+1
            comparative <- sceComparative[[clusterNTarget]]
            ## Comparative is in this case single cluster from sample_i+1
            textwide <- paste0("\r Comparing [cluster ", clusterNames[[sampleN]][clusterN], "] from dataset: ", sampleNames[sampleN], 
                               " with [cluster: ", clusterNames[[sampleNTarget]][clusterNTarget] , "] from dataset: ", sampleNames[sampleNTarget], spaces)
            message(textwide, appendLF=FALSE)
            ## Create a matrix A and B to compute the dot-product of all possible combinations of cluster comparison FoldChanges for each gene
            mat <- foldchangeComposition(root[features,], comparative[features,])
            ## We apply the mean, as the number of clusters between datasets could be different,
            ## just by doing the sum could be biased by the number of pair-wise FC comparisons
            matColmean <- colMeans(mat,na.rm=TRUE) # We obtain one value per gene
            topGenes <- head(colnames(mat)[order(matColmean, decreasing=TRUE)], n=topNFeatures)
            nNegative <- 0
            nPositive <- 0
            nNegative <- sum(matColmean< (-0.2), na.rm=TRUE)
            nPositive <- sum(matColmean> (0.20), na.rm=TRUE)
            sem <- round(sd(matColmean) / sqrt(nNegative + nPositive + 1), digits=2)
            ## Weight based on the number of concordant-discordant FCs:
            ## loosely based on cross-entropy
            weight <- round(log(max(abs(nPositive - nNegative), 1)), digits = 2) * sign(nPositive - nNegative) # max: Avoiding log 0
            similarity_weighted <- (sqrt(abs(sum(matColmean))) + weight) * sign(sum(matColmean) + weight)
            ## Save results
            for (gene in topGenes){
              results[nrow(results) + 1,] <- list(
                similarity_weighted, ## similarityValue
                sem, ## Standar Error of the Mean
                weight, ## Weight of the score
                sampleNames[sampleN], ## datasetL
                clusterNames[[sampleN]][clusterN], ## clusterL (left; source of comparison) -> clusterN corresponding to the loop
                sampleNames[sampleNTarget], ## datasetR
                clusterNames[[sampleNTarget]][clusterNTarget], ## clusterR (right; target of comparison) -> clusterNTarget corresponding to the internal loop
                gene, ## Top feature conserved
                matColmean[gene]) ## Feature Score
            }
          }
        if(is.infinite(topN)){
        results <- do.call("rbind", by(results, results$clusterL, function(x){x[order(x[,"similarityValue"], decreasing=TRUE),]}))
        }else{
        results <- do.call("rbind", by(results, results$clusterL, function(x){head(x[order(x[,"similarityValue"], decreasing=TRUE),], n=topN*topNFeatures)}))
        }
        summaryResults <- rbind.data.frame(summaryResults,results)
      }
    return(summaryResults)
    })
  }) # End main loop - sapply
  if(is(intermediateResults[[1]], "list")){
    # We obtain a list of list containing the data.frame results
    summaryResults <- do.call("rbind", lapply(intermediateResults, function(x){do.call(rbind, x)}))
  }else if(is(intermediateResults[[1]], "data.frame")){
    # Special case in which the number of cluster is the same on all datasets, we will get a data.frame
    summaryResults <- do.call("rbind", intermediateResults)
  }
  if(is.infinite(topN)){
  message("\n Ploting heatmap using the similarity values of clusters (topN=Inf).")
  show(similarityHeatmap(similarityTable=summaryResults))
  }else{
  message("\n Ploting graph using the similarity values of clusters.")
  plotClustersGraph(similarityTable=summaryResults)
  }
  message("Returning similarity table.")
  return(summaryResults)
}
