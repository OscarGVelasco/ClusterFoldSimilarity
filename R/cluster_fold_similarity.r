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
#' @param top_n Numeric. Specifies the number of target clusters with best similarity to report for each cluster comparison (default 1). If set to Inf, then all similarity values from all possible pairs of clusters are returned.
#' @param top_n_genes Numeric. Number of top genes that explains the clusters similarity to report for each cluster comparison (default 1). If top_n = Inf then top_n_genes is automatically set to 1.
#' 
#' @return The function returns a \linkS4class{DataFrame} containing the best top similarities between all possible pairs of single cell samples. Column values are:
#' \tabular{ll}{
#'    \code{similarity_value} \tab The top similarity value calculated between dataset_l:cluster_l and dataset_r. \cr
#'    \tab \cr
#'    \code{sem} \tab Standar Error of the Mean (SEM) of the mean of the values of the coeficient calculated for all genes. \cr
#'    \tab \cr
#'    \code{dataset_l} \tab Dataset left, the dataset/sample which has been used to be compared.  \cr
#'    \tab \cr
#'    \code{cluster_l} \tab Cluster left, the cluster source from dataset_l which has been compared. \cr
#'    \tab \cr
#'    \code{dataset_r} \tab Dataset right, the dataset/sample used for comparison against dataset_l. \cr
#'    \tab \cr
#'    \code{cluster_r} \tab Cluster right, the cluster target from dataset_r which is being compared with the cluster_l from dataset_l. \cr
#'    \tab \cr
#'    \code{top_gene_conserved} \tab The gene that showed most similar between cluster_l & cluster_r. \cr
#' }
#'
#' @examples 
#' library(Seurat)
#' library(SeuratData)
#' library(ClusterFoldSimilarity)
#' # load dataset
#' LoadData("ifnb")
#' # split the dataset into a list of two seurat objects (stim and CTRL)
#' singlecell.object.list <- SplitObject(ifnb, split.by = "stim")
#' # normalize, identify variable features and cluster for each dataset independently
#' singlecell.object.list <- lapply(X = singlecell.object.list, FUN = function(x){
#' x <- NormalizeData(x)
#' x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
#' x <- ScaleData(x,features = VariableFeatures(x))
#' x <- RunPCA(x, features = VariableFeatures(object = x))
#' x <- FindNeighbors(x, dims = 1:10)
#' x <- FindClusters(x, resolution = 0.1)
#' })
#' 
#' similarity.table <- cluster_fold_similarity(sce_list = singlecell.object.list)
#' 
#' @author Oscar Gonzalez-Velasco
#' @export
cluster_fold_similarity <- function(sce_list = NULL,
                                  sample_names = NULL,
                                  top_n = 1,
                                  top_n_genes = 1){
  if(is.null(sce_list) | (length(sce_list)<2)){
    stop("At least two Single Cell Experiments are needed for cluster comparison.")
  }
  if(is.infinite(top_n)){
    # if top_n clusters to report is set to Inf, then we return ALL similarity values from ALL the possible combination of clusters
    message("Returning similarity values from ALL pairs of clusters")
    message("  *(Note: by using top_n=Inf the top_n_genes is set automatically to 1)")
    top_n_genes = 1
  }
  is_seurat <- FALSE
  is_sce <- FALSE
  # Function starts by loading dependencies
  if(is(sce_list[[1]], 'Seurat')){
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package \"Seurat\" needed for this function to work. Please install it.",
           call. = FALSE)}
    is_seurat <- TRUE
    is_sce <- FALSE
  }
  if(is(sce_list[[1]], 'SingleCellExperiment')){
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.",
           call. = FALSE)}
    if(("normcounts" %in% names(sce_list[[1]]@assays)) == FALSE){
      stop("\"normcounts\" not found on the SingleCellExperiment object assays. \nPlease, set the normalized counts matrix using SingleCellExperiment::normcounts(sce_object) <- normcounts ",
           call. = FALSE)
    }
    is_seurat <- FALSE
    is_sce <- TRUE
  }
  if (!is_seurat & !is_sce){
    stop("One or more objects in the input list is neither of class Seurat nor SingleDataExperiment.")
  }
  summary_results <- data.frame(similarity_value=integer(),
                                sem=integer(),
                                dataset_l=integer(), 
                                cluster_l=integer(),
                                dataset_r=integer(),
                                cluster_r=integer(),
                                top_gene_conserved=character(),
                                stringsAsFactors=FALSE)
  # Select common genes in all samples 
  features <- Reduce(intersect,lapply(sce_list,function(x){rownames(x)}))
  if(length(features) == 0){
    stop("No common features between datasets. Please, select a subset of common variable features among all datasets.")
  }
  if(length(features) < 50){
    warning("The number of common features among datasets is: ",length(features),". More than 50 common features among all datasets is recomended.")
  }
  # Control filtered level factors using droplevels
  for (i in length(sce_list)){
    if(is_sce){
      if(is.null(SingleCellExperiment::colLabels(sce_list[[i]]))){
        stop("No clusters specified on colLabels on sample ",i)
      }
      if(is.factor(SingleCellExperiment::colLabels(sce_list[[i]]))){
        SingleCellExperiment::colLabels(sce_list[[i]]) <- droplevels(SingleCellExperiment::colLabels(sce_list[[i]]))
      }else{
        SingleCellExperiment::colLabels(sce_list[[i]]) <- factor(SingleCellExperiment::colLabels(sce_list[[i]]))
      }
      }
    if(is_seurat){
      if(is.factor(Seurat::Idents(sce_list[[i]]))){
        Seurat::Idents(sce_list[[i]]) <- droplevels(Seurat::Idents(sce_list[[i]]))
      }else{
        Seurat::Idents(sce_list[[i]]) <- factor(Seurat::Idents(sce_list[[i]]))
      }
      }
  }
  spaces <- "                  " # Empty spaces to clear the append of message report
  # Save the cluster name identification given by the user
  if(is_sce){cluster_names <- lapply(sce_list,function(x)levels(SingleCellExperiment::colLabels(x)))}
  if(is_seurat){cluster_names <- lapply(sce_list,function(x)levels(Seurat::Idents(x)))}
  # Check for dataset names:
  if(is.null(sample_names)){
    sample_names <- seq(length(sce_list))
  }else{
    if(!(length(sample_names) == length(sce_list))){
      stop("Number of names given in sample_names does not match the length of the single-cell experiment list.")
    }
  }
  # Calculate cluster FoldChange pairwise values:
  markers_sce_list <- list()
  # markers_sce_list <- lapply(sce_list,function(x){findMarkers(x,pval.type="all")})
  message("Computing fold changes.")
  if(is_sce){
    markers_sce_list <- lapply(sce_list,function(x){
      pairwise_cluster_fold_change(x = as.matrix(SingleCellExperiment::normcounts(x)),clusters = SingleCellExperiment::colLabels(x))
      })
  }
  if(is_seurat){
    markers_sce_list <- lapply(sce_list,function(x){
      pairwise_cluster_fold_change(x = as.matrix(Seurat::GetAssayData(x, slot = "data")),clusters = Seurat::Idents(x))
    })
  }
  # Main loop - samples
  for (i in seq_len(length(markers_sce_list))){
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
                              dataset_l=integer(), 
                              cluster_l=integer(),
                              dataset_r=integer(),
                              cluster_r=integer(),
                              top_gene_conserved=character(),
                              stringsAsFactors=FALSE)
          for(n in seq_along(sce_comparative)){
            # The sample for comparing will be the sample_i+1
            comparative <- sce_comparative[[n]]
            # Comparative = single cluster from sample_i+1
            textwide <- paste0("\r Comparing [cluster ",cluster_names[[i]][j],"] from dataset: ",sample_names[i]," with [cluster: ",cluster_names[[k]][n] ,"] from dataset: ",sample_names[k],spaces)
            message(textwide, appendLF=FALSE)
            # Create a matrix A and B to compute the dotproduct of all possible combinations of cluster comparison FoldChanges for each gene
            mat <- foldchange_composition(root[features,],comparative[features,])
            # We apply the mean, as the number of clusters between datasets could be different,
            # just by doing the sum could be biased by the number of pair-wise FC comparisons
            mat_colmean <- colMeans(mat,na.rm = TRUE) # We obtain one value per gene
            top_genes <- head(colnames(mat)[order(mat_colmean,decreasing = TRUE)],n=top_n_genes)
            n_negative = 0
            n_positive = 0
            n_negative <- sum(mat_colmean < 0,na.rm = TRUE)
            n_positive <- sum(mat_colmean > 0,na.rm = TRUE)
            diff_neg_pos <- n_negative - n_positive
            if(diff_neg_pos == 0 ){
              # if same n of neg. and pos. OR all are 0s we add +2 to avoid Inf and multiply by 1 (log2(2))
              diff_neg_pos=2
            }else if(abs(diff_neg_pos) == 1 ){
              # Num of pos. or neg. is exactly +1 larger -> multiply by 1 (log2(2)) to avoid multiply by 0
              diff_neg_pos=2
            }
            sem <- round(sd(mat_colmean) / sqrt(n_negative + n_positive), digits = 4)
            # Weight based on the number of concordant-discordant FCs:
            # log2 allows us to set 0 when we have equal number of + and -, and when - number is large, the weight is also
            # large, allowing us to penalize negative numbers as well
            weight <- round(log2(abs(diff_neg_pos)), digits = 2) # minimum possible value of weight is = 1
            similarity_weighted <- sqrt(abs(sum(mat_colmean)) * weight) * sign(sum(mat_colmean))
            # Save results
            for (g in top_genes){
              results[nrow(results)+1,] <- list(
                similarity_weighted, # Similarity_value
                sem, # Standar Error of the Mean
                # i, # dataset_l
                sample_names[i], # dataset_l
                cluster_names[[i]][j], # Cluster_l (left; source of comparison) -> j corresponding to the loop
                # k, # dataset_r
                sample_names[k], # dataset_r
                cluster_names[[k]][n], # Cluster_r (right; target of comparison) -> n corresponding to the internal loop
                g) # Top gene conserved
            }
          }
        if(is.infinite(top_n)){
        results <- data.table::rbindlist(by(results,results$cluster_l,function(x){x[order(x[,"similarity_value"],decreasing = TRUE),]}))
        }else{
        results <- data.table::rbindlist(by(results,results$cluster_l,function(x){head(x[order(x[,"similarity_value"],decreasing = TRUE),], n=top_n*top_n_genes)}))
        }
        summary_results <- rbind.data.frame(summary_results,results)
      }
    }
  }
  if(is.infinite(top_n)){
  message("\n Ploting heatmap using the similarity values of clusters (top_n=Inf).")
  print(similarity_heatmap(similarity_table = summary_results))
  }else{
  message("\n Ploting graph using the similarity values of clusters.")
  plot_clusters_graph(similarity.table = summary_results)
  }
  message("Returning similarity table.")
  return(summary_results)
}
