#' Calculate the dot product between all the possible combinations of foldchanges from diferent clusters.
#'
#' `foldchangeComposition()` returns a dataframe containing the best top similarities between all possible pairs of single cell samples.
#'
#' This function will perform the dot product of each possible combination of foldchanges, by constructing two dataframes: 
#' one with the source cluster's foldchanges and the other with the foldchange values of a target sample's cluster.
#' The computation of all the possible combinations is the hadamard product of the matrix.
#'
#' @param root Dataframe. Foldchanges between a source cluster and all the other clusters found on a sample.
#' @param comparative Dataframe. Foldchanges between a cluster and all the other clusters found on a second sample to be compared with the (root) cluster foldchanges.
#' 
#' @return A dataframe containing the hadamard product of all the possible combinations of foldchanges.
#' 
#' @keywords internal
foldchangeComposition <- function(root=NULL, comparative=NULL){
  ## Avoiding 1 column matrix to vector conversion when sub-setting index:
  root <- as.matrix(root)
  comparative <- as.matrix(comparative)
  tmp <- t(matrix(unlist(root), nrow=nrow(root), ncol=ncol(root)*ncol(comparative))) * matrix(rep(t(comparative), each=ncol(root)), ncol=nrow(root))
  colnames(tmp) <- rownames(root)
  return(tmp)
  }
