#' Calculate the gene mean expression Fold Change between all possible combinations of clusters.
#'
#' `pairwiseClusterFoldChange()` returns a list of dataframes containing the pairwise fold changes between all combinations of cluster.
#'
#' This function will perform fold change estimation from the mean feature´s expression between all possible combination of clusters specified on colLabels inside the sc object.
#' Bayesian Estimation of FoldChanges and Pseudocounts adapted from:
#' Florian Erhard, Estimating pseudocounts and fold changes for digital expression measurements, Bioinformatics, 
#'   Volume 34, Issue 23, December 2018, Pages 4054–4063, https://doi.org/10.1093/bioinformatics/bty471
#' Please consider citing also Erhard et. al. paper when using ClusterFoldChange.
#'
#' @param x Dataframe. Normalized counts containint gene expression.
#' @param clusters Factor. A vector of corresponding cluster for each sample of (x).
#' @param nSubsampling Numeric. Number of random samplings of cells to achieve fold change stability.
#' 
#' @return A list of dataframes containing the pairwise fold changes between all combinations of cluster.
#' 
#' @author Oscar Gonzalez-Velasco
#' @importFrom("methods", "is")
#' @importFrom("stats", "median", "optim", "pnorm", "quantile", "var")
#' @keywords internal
pairwiseClusterFoldChange <- function(x, clusters, nSubsampling){
  countData <- x
  minNumFeatures <- 3
  cellPortion <- 1/3
  fcFunct <- function(x,y) "-"(log2(x),log2(y)) 
  clusters <- as.factor(clusters)
  cellGroups <- split(seq_len(ncol(countData)), clusters)
  j <- expand.grid(levels(clusters), levels(clusters))
  j <- t(j[!(j[,1] == j[,2]),])[c(2, 1),]
  # Sub-sampling of cells - we will randomly select "nSubsampling" times the "cellPortion" of the cells and calculate the fold change mean
  listOfFolds <- lapply(seq(1, nSubsampling), function(sampling){
    i <- lapply(cellGroups,function(cells)sample(x=cells, size=max((length(cells)*cellPortion), 1)))
    x <- vapply(i, function(i){ rowMeans(countData[,i,drop=FALSE])}, FUN.VALUE=double(nrow(countData)))
    logFolds <- mapply(function(i, k){
        numerator <- x[,i]
        denominator <- x[,k]
        informativeFeatures <- (numerator > 0) & (denominator > 0) # Filter only genes with exprs. in BOTH conditions as informative
        if(sum(informativeFeatures) >= minNumFeatures){ # A minimum of 3 features are needed
          numerator <- numerator[informativeFeatures]
          denominator <- denominator[informativeFeatures]
          folds <- fcFunct(numerator, denominator)
          muPri <- mean(folds) # We choose mean as prior
          varPri <- max((quantile(folds, pnorm(1))-muPri)^2, (-quantile(folds, pnorm(-1))+muPri)^2) # pnorm(i) probability in a normal distr. to find the FC i (mean of FC dist. = 0)
        }
      else{
        muPri <- Inf # If not enough informative features we apply the pseudocount approach
        varPri <- Inf
      }
      if (is.infinite(muPri) || is.infinite(varPri)) {
        # pseudocount = 1 approach
        muPri <- mean(fcFunct(x[,i]+1, x[,k]+1))
        varPri <- var(fcFunct(x[,i]+1, x[,k]+1))
      }
      sqrFunct <- function(params)((digamma(params[1])-digamma(params[2])-muPri)^2 + (trigamma(params[1])+trigamma(params[2])-varPri)^2)
      pseudoCounts <- optim(c(1, 1), sqrFunct)$par
      finalfolds <- (digamma(x[,i]+pseudoCounts[1]) - digamma(x[,k]+pseudoCounts[2])) / log(2)
      return(finalfolds)
      }, j[1,], j[2,])
  })
  # Calculate the element-wise mean of the sub-samplings of fold changes
  finalLogFolds <- Reduce(`+`, listOfFolds) / length(listOfFolds)
  # Set the extremely small fold changes to 0
  #finalLogFolds[finalLogFolds<0.1 & finalLogFolds>-0.1] <- 0
  finalLogFolds <- apply(finalLogFolds, 2,function(z)(z - median(z))) # Normalization by median
  colnames(finalLogFolds) <- paste("logFC", j[1,], j[2,], sep='.')
  n <- length(levels(clusters)) - 1
  nr <- ncol(finalLogFolds)
  fcValues <- lapply(split.data.frame(t(finalLogFolds), rep(seq_len(ceiling(nr/n)), each=n, length.out=nr)), t) # A list of data.frames per cluster
  return(fcValues)
}
