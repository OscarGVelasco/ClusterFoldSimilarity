#' Calculate cluster similarity between clusters from different single cell samples.
#'
#' `clusterFoldSimilarity()` returns a dataframe containing the best top similarities between all possible pairs of single cell samples.
#'
#' This function will calculate a similarity coeficient using the fold changes of shared features (e.g.:genes for a single-cell RNA-Seq, peaks for ATAC-Seq) among clusters, or user-defined-groups, from different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster/group i from sample n and all the target clusters/groups in sample j.
#'
#' @param sceList List. A list of Single Cell Experiments or Seurat objects. At least 2 are needed. The objects are expected to have cluster or label groups set as identity class.
#' @param sampleNames Character Vector. Specify the sample names, if not a number corresponding with its position on (sceList).
#' @param topN Numeric. Specifies the number of target clusters with best similarity to report for each cluster comparison (default 1). If set to Inf, then all similarity values from all possible pairs of clusters are returned.
#' @param topNFeatures Numeric. Number of top features that explains the clusters similarity to report for each cluster comparison (default 1). If topN = Inf then topNFeatures is automatically set to 1.
#' @param nSubsampling Numeric. Number of random sampling of cells to achieve fold change stability (default 15).
#' @param parallel Boolean. Whether to use parallel computing using BiocParallel or not (default FALSE).
#' 
#' @return The function returns a DataFrame containing the best top similarities between all possible pairs of single cell samples. Column values are:
#' \tabular{ll}{
#'    \code{similarityValue} \tab The top similarity value calculated between datasetL:clusterL and datasetR. \cr
#'    \tab \cr
#'    \code{sem} \tab Standar Error of the Mean (SEM) of the mean of the values of the coeficient calculated for all genes. \cr
#'    \tab \cr
#'    \code{w} \tab Weight associated with the score value. \cr
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
#' }
#'
#' @examples 
#' if (requireNamespace("Seurat") & requireNamespace("SingleCellExperiment") & requireNamespace("scRNAseq") & requireNamespace("SeuratObject")){
#' library(ClusterFoldSimilarity)
#' # Mouse brain single-cell RNA-seq 1 from Romanov et. al.
#' mouse.brain.romanov <- scRNAseq::RomanovBrainData(ensembl = TRUE)
#' colnames(mouse.brain.romanov) <- colData(mouse.brain.romanov)$cellID
#' rownames(colData(mouse.brain.romanov)) <- colData(mouse.brain.romanov)$cellID
#' singlecell.1.seurat <- Seurat::CreateSeuratObject(counts = counts(mouse.brain.romanov),meta.data = as.data.frame(colData(mouse.brain.romanov)))
#' 
#' # Mouse brain single-cell RNA-seq 2 from Zeisel et. al.
#' mouse.brain.zei <- scRNAseq::ZeiselBrainData(ensembl = TRUE)
#' singlecell.2.seurat <- Seurat::CreateSeuratObject(counts = counts(mouse.brain.zei),meta.data = as.data.frame(colData(mouse.brain.zei)))
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
#' # Calculate the similarity between the two datasets:
#' similarity.table <- clusterFoldSimilarity(sceList=singlecell.object.list, topN=1, nSubsampling=18)
#' head(similarity.table)
#' }
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom BiocParallel bplapply
#' @export
clusterFoldSimilarity <- function(sceList=NULL, sampleNames=NULL, topN=1, topNFeatures=1, nSubsampling=15, parallel=FALSE){
  if(is.null(sceList) | (length(sceList)<2)){
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
  if(is(sceList[[1]], 'Seurat')){
    if (!requireNamespace("Seurat", quietly=TRUE)) {
      stop("Package \"Seurat\" needed for this function to work. Please install it.",
           call.=FALSE)}
    isSeurat <- TRUE
    isSce <- FALSE
  }
  if(is(sceList[[1]], 'SingleCellExperiment')){
    if (!requireNamespace("SingleCellExperiment", quietly=TRUE)) {
      stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.",
           call.=FALSE)}
    # if(("normcounts" %in% names(sceList[[1]]@assays)) == FALSE){
    #   stop("\"normcounts\" not found on the SingleCellExperiment object assays. \nPlease, set the normalized counts matrix using SingleCellExperiment::normcounts(sce_object) <- normcounts ",
    #        call.=FALSE)
    # }
    isSeurat <- FALSE
    isSce <- TRUE
  }
  if (!isSeurat & !isSce){
    stop("One or more objects in the input list is neither of class Seurat nor SingleDataExperiment.")
  }
  functToApply <- base::lapply
  if(isTRUE(parallel)){
    #BiocParallel::MulticoreParam(workers = 6)
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
                                geneScore=integer(),
                                stringsAsFactors=FALSE)
  ## Select common genes in all samples 
  features <- Reduce(intersect,lapply(sceList,function(x){rownames(x)}))
  if(length(features) == 0){
    stop("No common features between datasets. Please, select a subset of common variable features among all datasets.")
  }
  if(length(features) < 50){
    warning("The number of common features among datasets is: ",length(features),". More than 50 common features among all datasets is recomended.")
  }
  ## Control filtered level factors using droplevels
  for (i in length(sceList)){
    if(isSce){
      if(is.null(SingleCellExperiment::colLabels(sceList[[i]]))){
        stop("No clusters specified on colLabels on sample ",i)
      }
      if(is.factor(SingleCellExperiment::colLabels(sceList[[i]]))){
        SingleCellExperiment::colLabels(sceList[[i]]) <- droplevels(SingleCellExperiment::colLabels(sceList[[i]]))
      }else{
        SingleCellExperiment::colLabels(sceList[[i]]) <- factor(SingleCellExperiment::colLabels(sceList[[i]]))
      }
      }
    if(isSeurat){
      if(is.factor(Seurat::Idents(sceList[[i]]))){
        Seurat::Idents(sceList[[i]]) <- droplevels(Seurat::Idents(sceList[[i]]))
      }else{
        Seurat::Idents(sceList[[i]]) <- factor(Seurat::Idents(sceList[[i]]))
      }
      }
  }
  spaces <- "                  " ## Empty spaces to clear the append of message report
  ## Save the cluster name identification given by the user
  if(isSce){clusterNames <- lapply(sceList,function(x)levels(SingleCellExperiment::colLabels(x)))}
  if(isSeurat){clusterNames <- lapply(sceList,function(x)levels(Seurat::Idents(x)))}
  ## Check for dataset names:
  if(is.null(sampleNames)){
    sampleNames <- seq(length(sceList))
  }else{
    if(!(length(sampleNames) == length(sceList))){
      stop("Number of names given in sampleNames does not match the length of the single-cell experiment list.")
    }
  }
  ## Calculate cluster FoldChange pairwise values:
  markersSceList <- list()
  message(paste("Using a common set of", length(features), "features."))
  # Testing the number of draws/subsamplings of cells needed to see *all cells*
  cellDraw <- function(n){
    percentage = (1/3) * n
    return(sum(rep(1, n+1) / seq(1, n+1)) * (n/percentage))
  }
  if(isSce){
    nOfDraws <- unlist(sapply(1:length(sceList),function(i)sapply(c(table(SingleCellExperiment::colLabels(sceList[[i]]))),cellDraw)))
  }
  if(isSeurat){
    nOfDraws <- unlist(sapply(1:length(sceList),function(i)sapply(c(table(Seurat::Idents(sceList[[i]]))),cellDraw)))
  }
  # We select as optimal the subsampling n percentile 85% across all cell groups
  message(paste0("Using a cell subsampling of n=",nSubsampling," (recomended n=",round(quantile(nOfDraws, probs=0.85)),")"))
  message("Computing fold changes.")
  if(isSce){
    markersSceList <- functToApply(sceList, function(x){
      pairwiseClusterFoldChange(x=as.matrix(SingleCellExperiment::counts(x)), clusters=SingleCellExperiment::colLabels(x), nSubsampling=nSubsampling, functToApply=functToApply) # Raw counts
    })
  }
  if(isSeurat){
    markersSceList <- functToApply(sceList, function(x){
      pairwiseClusterFoldChange(x=as.matrix(Seurat::GetAssayData(x, slot="counts")), clusters=Seurat::Idents(x), nSubsampling=nSubsampling, functToApply=functToApply) # Raw counts from Seurat object
    })
  }
  ## Main loop - samples
  intermediateResults <- sapply(seq_len(length(markersSceList)),function(i){
    ## Pick a sample in ascending order
    y <- markersSceList[[i]]
    intermediateResults <- functToApply(seq_along(y),function(j){ # Parallel computing if selected
      ## We choose a root cluster and compare it with clusters of samples i+1+...+n
      root <- y[[j]]
      for (k in seq(from=1, to=length(markersSceList), by=1)){
        if(k == i)next() # If root and target samples are the same, we go for the next sample and skip this loop
        sceComparative <- markersSceList[[k]]
        results <- data.frame(similarityValue=integer(),
                              sem=integer(),
                              w=integer(),
                              datasetL=integer(), 
                              clusterL=integer(),
                              datasetR=integer(),
                              clusterR=integer(),
                              topFeatureConserved=character(),
                              geneScore=integer(),
                              stringsAsFactors=FALSE)
          for(n in seq_along(sceComparative)){
            ## The sample for comparing will be the sample_i+1
            comparative <- sceComparative[[n]]
            ## Comparative is in this case single cluster from sample_i+1
            textwide <- paste0("\r Comparing [cluster ",clusterNames[[i]][j],"] from dataset: ",sampleNames[i]," with [cluster: ",clusterNames[[k]][n] ,"] from dataset: ",sampleNames[k],spaces)
            message(textwide, appendLF=FALSE)
            ## Create a matrix A and B to compute the dotproduct of all possible combinations of cluster comparison FoldChanges for each gene
            mat <- foldchangeComposition(root[features,],comparative[features,])
            ## We apply the mean, as the number of clusters between datasets could be different,
            ## just by doing the sum could be biased by the number of pair-wise FC comparisons
            matColmean <- colMeans(mat,na.rm=TRUE) # We obtain one value per gene
            topGenes <- head(colnames(mat)[order(matColmean, decreasing=TRUE)],n=topNFeatures)
            nNegative <- 0
            nPositive <- 0
            nNegative <- sum(matColmean< (-0.2), na.rm=TRUE)
            nPositive <- sum(matColmean> (0.20), na.rm=TRUE)
            sem <- round(sd(matColmean) / sqrt(nNegative + nPositive), digits=2)
            ## Weight based on the number of concordant-discordant FCs:
            ## loosely based on cross-entropy
            weight <- round(log(max(abs(nPositive - nNegative), 1)), digits = 2) * sign(nPositive - nNegative) # max: Avoiding log 0
            similarity_weighted <- (sqrt(abs(sum(matColmean))) + weight) * sign(sum(matColmean) + weight)
            ## Save results
            for (g in topGenes){
              results[nrow(results) + 1,] <- list(
                similarity_weighted, ## similarityValue
                sem, ## Standar Error of the Mean
                weight, ## Weight of the score
                sampleNames[i], ## datasetL
                clusterNames[[i]][j], ## clusterL (left; source of comparison) -> j corresponding to the loop
                sampleNames[k], ## datasetR
                clusterNames[[k]][n], ## clusterR (right; target of comparison) -> n corresponding to the internal loop
                g, ## Top feature conserved
                matColmean[g]) ## Gene Score
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
  summaryResults <- do.call("rbind",lapply(intermediateResults,function(x)do.call(rbind, x)))
  if(is.infinite(topN)){
  message("\n Ploting heatmap using the similarity values of clusters (topN=Inf).")
  print(similarityHeatmap(similarityTable=summaryResults))
  }else{
  message("\n Ploting graph using the similarity values of clusters.")
  plotClustersGraph(similarityTable=summaryResults)
  }
  message("Returning similarity table.")
  return(summaryResults)
}
