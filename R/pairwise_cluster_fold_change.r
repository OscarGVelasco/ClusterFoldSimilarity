#' Calculate the gene mean expression Fold Change between all possible combinations of clusters.
#'
#' `pairwiseClusterFoldChange()` returns a list of dataframes containing the pairwise fold changes between all combinations of cluster.
#'
#' This function will perform fold change estimation from the mean feature´s expression between all possible combination of clusters specified on colLabels inside the sc object.
#' Bayesian Estimation of FoldChanges and Pseudocounts adapted from:
#' Florian Erhard, Estimating pseudocounts and fold changes for digital expression measurements, Bioinformatics, 
#'   Volume 34, Issue 23, December 2018, Pages 4054–4063, https://doi.org/10.1093/bioinformatics/bty471
#' Please consider citing also Erhard et. al. paper when using ClusterFoldSimilarity.
#'
#' @param countData Matrix. Normalized counts containing gene expression.
#' @param clusters Factor. A vector of corresponding cluster for each sample of (x).
#' @param nSubsampling Numeric. Number of random sampling of cells to achieve fold change stability.
#' 
#' @return A list of dataframes containing the pairwise fold changes between all combinations of cluster.
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom methods is
#' @importFrom stats median optim pnorm quantile var
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix rowMeans
#' @keywords internal
pairwiseClusterFoldChange <- function(countData, clusters, nSubsampling, functToApply){
  minNumFeatures <- 3
  cellPortion <- 1/3
  fcFunct <- function(x, y) "-"(log2(x), log2(y))
  # To avoid integer overflow use counts^2 instead of counts*counts:
  pseudoCount <- function(counts){counts + sqrt((counts^2)+1)}
  clusters <- as.factor(clusters)
  cellGroups <- split(seq_len(ncol(countData)), clusters)
  groupsByCluster <- expand.grid(levels(clusters), levels(clusters))
  groupsByCluster <- t(groupsByCluster[!(groupsByCluster[,1] == groupsByCluster[,2]),])[c(2, 1),]
  # Sub-sampling of cells - we will randomly select "nSubsampling" times the "cellPortion" of the cells and calculate the fold change mean
  listOfFolds <- functToApply(seq(1, nSubsampling), function(sampling){
    i <- lapply(cellGroups, function(cells)sample(x=cells, size=max((length(cells)*cellPortion), 1)))
    meanCounts <- vapply(i, function(i){ Matrix::rowMeans(countData[,i,drop=FALSE])}, FUN.VALUE=double(nrow(countData)))
    logFolds <- mapply(function(group1, group2){
        numerator <- meanCounts[, group1]
        denominator <- meanCounts[, group2]
        informativeFeatures <- (numerator > 0) & (denominator > 0) # Filter only genes with exprs. in BOTH conditions as informative
        if(sum(informativeFeatures) >= minNumFeatures){ # A minimum of 3 features are needed
          numerator <- numerator[informativeFeatures]
          denominator <- denominator[informativeFeatures]
          folds <- fcFunct(numerator, denominator)
          muPri <- mean(folds) # We choose mean as prior
          varPri <- max((quantile(folds, pnorm(1))-muPri)^2, (-quantile(folds, pnorm(-1))+muPri)^2) # pnorm(group1) probability in a normal distr. to find the FC group1 (mean of FC dist. = 0)
        }
      else{
        muPri <- Inf # If not enough informative features we apply the pseudocount approach
        varPri <- Inf
      }
      if (is.infinite(muPri) || is.infinite(varPri)) {
        muPri <- mean(fcFunct(pseudoCount(meanCounts[, group1]), pseudoCount(meanCounts[, group2])))
        varPri <- var(fcFunct(pseudoCount(meanCounts[, group1]), pseudoCount(meanCounts[, group2])))
      }
      sqrFunct <- function(params)((digamma(params[1])-digamma(params[2])-muPri)^2 + (trigamma(params[1])+trigamma(params[2])-varPri)^2)
      pseudoCounts <- optim(c(1, 1), sqrFunct)$par
      finalfolds <- (digamma(meanCounts[, group1]+pseudoCounts[1]) - digamma(meanCounts[, group2]+pseudoCounts[2])) / log(2)
      return(finalfolds)
      }, groupsByCluster[1,], groupsByCluster[2,])
  })
  # Calculate the element-wise mean of the sub-samplings of fold changes
  finalLogFolds <- Reduce(`+`, listOfFolds) / length(listOfFolds)
  finalLogFolds <- apply(finalLogFolds, 2,function(z)(z - stats::median(z))) # Normalization by median
  colnames(finalLogFolds) <- paste("logFC", groupsByCluster[1,], groupsByCluster[2,], sep='.')
  n <- length(levels(clusters)) - 1
  nr <- ncol(finalLogFolds)
  fcValues <- lapply(split.data.frame(t(finalLogFolds), rep(seq_len(ceiling(nr/n)), each=n, length.out=nr)), t) # A list of data.frames per cluster
  return(fcValues)
}
