#' Calculate the gene mean expression Fold Change between all possible combinations of clusters.
#'
#' `pairwise_cluster_fold_change()` returns a list of dataframes containing the pairwise fold changes between all combinations of cluster.
#'
#' This function will perform fold change ratios of mean gene expression between all possible combination of clusters specified on colLabels inside the sc object.
#'
#' @param x Dataframe. Normalized counts containint gene expression.
#' @param comparative Factor. A vector of corresponding cluster for each sample of (x).
#' 
#' @return A list of dataframes containing the pairwise fold changes between all combinations of cluster.
#' 
#' @author Oscar Gonzalez-Velasco
#' @keywords internal
pairwise_cluster_fold_change <- function(x, clusters){
  divis_funct <- function(x,y) "/"(x,y) 
  clusters <- as.factor(clusters)
  i <- split(seq_len(ncol(x)), clusters)
  x <- vapply(i, function(i){ rowMeans(x[,i])},FUN.VALUE = double(nrow(x)))
  j <- expand.grid(levels(clusters),levels(clusters))
  j <- t(j[!(j[,1] == j[,2]),])[c(2,1),]
  log.folds <- log2(divis_funct(x[,j[1,]]+1, x[,j[2,]]+1))
  colnames(log.folds) <- paste("logFC",j[1,], j[2,], sep = '.')
  n <- length(levels(clusters))-1
  nr <- ncol(log.folds)
  fc.values <- lapply(split.data.frame(t(log.folds), rep(seq_len(ceiling(nr/n)), each=n, length.out=nr)),t)
  return(fc.values)
}
