#' Calculate cluster similarity between clusters from different single cell samples.
#'
#' `cluster_fold_similarity()` returns a dataframe containing the best top similarities between all possible pairs of single cell samples.
#'
#' This function will calculate a similarity coeficient using the fold changes of shared genes among clusters of different samples/batches. The similarity coeficient
#' is calculated using the dotproduct of every pairwise combination of Fold Changes between a source cluster i of sample n and all the target clusters in sample j.
#'
#' @param sce_list List. A list of single cell experiments. At least 2 single cell experiments are needed. The objects are expected to be subsets of the original data, 
#' containing the selected features (e.g. 2000 most variable genes).
#' @param sample_names Character Vector. Specify the sample names, if not a number corresponding with its position on (sce_list).
#' @param top_n Numeric. Specifies the number of target clusters with best similarity to report for each cluster comparison (default 1).
#' @param top_n_genes Numeric. Number of top genes that explains the clusters similarity to report for each cluster comparison (default 1).
#' 
#' @return The function returns a \linkS4class{DataFrame} containing the best top similarities between all possible pairs of single cell samples. Column values are:
#' \tabular{ll}{
#'    \code{similarity_value} \tab The top similarity value calculated between sample_l:cluster_l and sample_r. \cr
#'    \tab \cr
#'    \code{sem} \tab Standar Error of the Mean (SEM) of the mean of the values of the coeficient calculated for all genes. \cr
#'    \tab \cr
#'    \code{sample_l} \tab Sample left, the sample source which is being compared.  \cr
#'    \tab \cr
#'    \code{cluster_l} \tab Cluster left, the cluster source which is being compared. \cr
#'    \tab \cr
#'    \code{sample_r} \tab Sample right, the sample target which is being compared with sample_l. \cr
#'    \tab \cr
#'    \code{cluster_r} \tab Cluster right, the cluster target which is being compared with cluster_l. \cr
#'    \tab \cr
#'    \code{top_gene_conserved} \tab The gene that showed most similar between cluster_l & cluster_r. \cr
#' }
#' 
#' @export
cluster_fold_similarity <- function(sce_list = NULL,
                                  sample_names = NULL,
                                  top_n = 1,
                                  top_n_genes = 1){
  if(is.null(sce_list) | (length(sce_list)<2)){
    stop("At least two Single Cell Experiments are needed for cluster comparison.")
  }
  is_seurat <- FALSE
  is_sce <- FALSE
  # Function starts by loading dependencies
  if(class(sce_list[[1]]) == "Seurat" ){
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package \"Seurat\" needed for this function to work. Please install it.",
           call. = FALSE)}
    is_seurat <- TRUE
    is_sce <- FALSE
  }
  if(class(sce_list[[1]]) == "SingleCellExperiment" ){
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.",
           call. = FALSE)}
    is_seurat <- FALSE
    is_sce <- TRUE
  }
  if (!is_seurat & !is_sce){
    stop("One or more objects in the input list is neither of class Seurat nor SingleDataExperiment.")
  }
  summary_results <- data.frame(similarity_value=integer(),
                                sem=integer(),
                                sample_l=integer(), 
                                cluster_l=integer(),
                                sample_r=integer(),
                                cluster_r=integer(),
                                top_gene_conserved=character(),
                                stringsAsFactors=FALSE)
  # Select common genes in all samples 
  features <- Reduce(intersect,lapply(sce_list,function(x){rownames(x)}))
  if(length(features) == 0){
    stop("No common genes between samples. Please, select a subset of common variable genes among all samples.")
  }
  if(length(features) < 100){
    warning(paste("The number of common genes among samples is: ",length(features),". More than 100 common genes among all samples is recomended.",sep = ""))
  }
  # Control filtered level factors using droplevels
  for (i in length(sce_list)){
    if(is_sce){
      if(is.null(SingleCellExperiment::colLabels(sce_list[[i]]))){
        stop(paste("No clusters specified on colLabels on sample",i))
      }
      SingleCellExperiment::colLabels(sce_list[[i]]) <- droplevels(SingleCellExperiment::colLabels(sce_list[[i]]))
    }
    if(is_seurat){
      Seurat::Idents(sce_list[[i]]) <- droplevels(Seurat::Idents(sce_list[[i]]))
    }
  }
  # Save the cluster name identification given by the user
  if(is_sce){cluster_names <- lapply(sce_list,function(x)levels(SingleCellExperiment::colLabels(x)))}
  if(is_seurat){cluster_names <- lapply(sce_list,function(x)levels(Seurat::Idents(x)))}
  # Calculate cluster FoldChange pairwise values:
  markers_sce_list <- list()
  # markers_sce_list <- lapply(sce_list,function(x){findMarkers(x,pval.type="all")})
  message("Computing fold changes.")
  if(is_sce){
    markers_sce_list <- lapply(sce_list,function(x){
      pairwise_cluster_fold_change(x = as.matrix(SingleCellExperiment::logcounts(x)),clusters = SingleCellExperiment::colLabels(x))
      })
  }
  if(is_seurat){
    markers_sce_list <- lapply(sce_list,function(x){
      pairwise_cluster_fold_change(x = as.matrix(Seurat::GetAssayData(x, slot = "data")),clusters = Seurat::Idents(x))
    })
  }
  # Main loop - samples
  for (i in 1:(length(markers_sce_list))){
    # Pick a sample in ascending order
    y <- markers_sce_list[[i]]
    for (j in seq_along(y)){
      # We choose a root cluster and compare it with clusters of samples i+1+...+n
      root <- y[[j]]
      # for (k in seq(from=i+1,to=length(markers_sce_list),by = 1)){
      for (k in seq(from=1,to=length(markers_sce_list),by = 1)){
        if(k == i)next() # If root and target samples are the same, we go for the next sample and skip this loop
        sce_comparative <- markers_sce_list[[k]]
        results <- data.frame(similarity_value=integer(),
                              sem=integer(),
                              sample_l=integer(), 
                              cluster_l=integer(),
                              sample_r=integer(),
                              cluster_r=integer(),
                              top_gene_conserved=character(),
                              stringsAsFactors=FALSE)
          for(n in seq_along(sce_comparative)){
            # The sample for comparing will be the sample_i+1
            comparative <- sce_comparative[[n]]
            # Comparative = single cluster from sample_i+1
            message(paste("Comparing [ cluster",cluster_names[[i]][j],"] from sample:",i,"with [ cluster:",cluster_names[[k]][n] ,"] from sample:",k))
            # Create a matrix A and B to compute the dotproduct of all possible combinations of cluster comparison FoldChanges for each gene
            mat <- foldchange_composition(root[features,],comparative[features,])
            mat_colmean <- colMeans(mat)
            top_genes <- head(colnames(mat)[order(mat_colmean,decreasing = T)],n=top_n_genes)
            n_negative <- sum(mat_colmean < 0)
            n_positive <- sum(mat_colmean > 0)
            sem <- round(sd(mat_colmean) / sqrt(n_negative + n_positive), digits = 4)
            weight <- round( sqrt(abs(n_negative - n_positive) / (n_negative + n_positive)), digits = 4) # Weight based on the number of concordant-discordant FCs
            # We divide the similarity by the SEM value -> then we scale it using sqrt and conserve the sign
            # similarity_weighted <- sqrt(abs(sum(mat_colmean)) / sem) * sign(sum(mat_colmean))
            similarity_weighted <- sqrt(abs(sum(mat_colmean))) * weight * sign(sum(mat_colmean))
            # Save results
            for (g in top_genes){
              results[nrow(results)+1,] <- list(
                similarity_weighted, # Similarity_value
                sem, # Standar Error of the Mean
                i, # Sample_l
                cluster_names[[i]][j], # Cluster_l (left; source of comparison) -> j corresponding to the loop
                k, # Sample_r
                cluster_names[[k]][n], # Cluster_r (right; target of comparison) -> n corresponding to the internal loop
                g) # Top gene conserved
            }
          }
        results <- data.table::rbindlist(by(results,results$cluster_l,function(x){head(x[order(x[,"similarity_value"],decreasing = T),], n=top_n*top_n_genes)}))
        summary_results <- rbind.data.frame(summary_results,results)
      }
    }
  }
  message("Returning similarity table.")
  return(summary_results)
}
