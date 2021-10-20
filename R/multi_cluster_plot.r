#' Plot several single-cell datasets/samples dimensionality reductions (PCA, UMAP, tSNE) in a single plot 3D plot, each layer corresponding to a dataset/sample.
#'
#' `multi_cluster_plot()` plot dimensionality reduction (PCA, UMAP, tSNE) of a batch of samples/datasets.
#'
#' This function will plot a list of single-cell datasets/samples dimensionality reduction in a 3D layered plot. 
#'
#' @param sce_list List. A list of single cell experiments. At least 2 single cell experiments are needed.
#' @param dim_method Character. Specify the dimensionality reduction used, either PCA, UMAP or tSNE. Default: PCA.
#' @param pdf_name Character. If specify the plot is saved on a PDF file with the given name. Default: plot in the current graphic device.
#' 
#' @return A dimensionality reduction 3D plot (PCA, UMAP, tSNE) of a batch of samples/datasets.
#' 
#' @import scatterplot3d
#' 
#' @export
multi_cluster_plot <- function(sce_list = NULL,dim_method="PCA",pdf_name = NULL){

  if(is.null(sce_list) | (length(sce_list)<2)){
    stop("At least two Single Cell Experiments are needed for cluster comparison.")
  }
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
  if(dim_method == "PCA"){
    dimMeth <- "PCA"
    if(is_sce){useDimFunc <- scater::calculatePCA}
    if(is_seurat){useDimFunc <- Seurat::RunPCA}
  }
  if(dim_method == "UMAP"){
    dimMeth <- "UMAP"
    if(is_sce){useDimFunc <- scater::calculateUMAP}
    if(is_seurat){useDimFunc <- Seurat::RunUMAP}
      }
  if(dim_method == "tSNE"){
    dimMeth <- "tSNE"
    if(is_sce){useDimFunc <- scater::calculateTSNE}
    if(is_seurat){useDimFunc <- Seurat::RunTSNE}
  }
  
  if(is_seurat){getDimFunc <- function(x,dim_method = dimMeth,useDimFunc=useDimFunc){
    if(dimMeth %in% names(x@reductions)){return(Embeddings(x[[dim_method]])[,1:2]) # If the dimensionality has been previously calculated, we obtain the values.
    }else{
        return(useDimFunc(x)) # If the required dimensionality is NOT in the object, we calculate it.
      }}
  }
  if(is_sce){getDimFunc <- function(x,dim_method = dimMeth,useDimFunc=useDimFunc){
    if(dimMeth %in% reducedDimNames(x)){return(reducedDim(x,dimMeth)[,1:2]) # If the dimensionality has been previously calculated, we obtain the values.
    }else{
      return(useDimFunc(x)) # If the required dimensionality is NOT in the object, we calculate it.
    }}
  }
  
  # if(is_sce){
  #   if(dimMeth %in% reducedDimNames(sce_list[[1]])){data <- reducedDim(sce_list[[1]],dimMeth)}else{data <- useDimFunc(sce_list[[1]])}
  #   }
  
  data <- cbind.data.frame(getDimFunc(sce_list[[1]]),sample_id=1, cluster = colLabels(sce_list[[1]]))

    for (i in seq_along(sce_list)[-1]){
    data <- rbind.data.frame(data,cbind.data.frame(getDimFunc(sce_list[[i]]),sample_id=i, cluster = colLabels(sce_list[[i]])))
    # if(dimMeth %in% reducedDimNames(sce_list[[i]])){
    #   data <- rbind.data.frame(data,cbind.data.frame(reducedDim(sce_list[[i]],dimMeth),sample_id=i, cluster = colLabels(sce_list[[i]])))
    #  }else{
    #   data <- rbind.data.frame(data,cbind.data.frame(useDimFunc(sce_list[[i]]),sample_id=i, cluster = colLabels(sce_list[[i]])))
    #   }
  }
  
  data$sample_id <- data$sample_id
  y_max <- max(data[,1]) +2
  y_min <- min(data[,1]) -2
  z_max <- max(data[,2]) +2
  z_min <- min(data[,2]) -2
  if(!is.null(pdf_name)){
    graphics.off()
    pdf(file = pdf_name, width = 40, height = 16)
  }
  
  # the_colors <- ifelse((as.numeric(data$sample_id) %% 2) == 0,0.6,1)
  # adjustcolor(col = as.numeric(data$cluster),alpha.f = the_colors)
  spd <- scatterplot3d::scatterplot3d(x = data$sample_id,y = data[,1],z = data[,2],color = data$cluster,
                                      angle = 60,pch = 19,box = FALSE,col.grid = "grey",
                                      cex.symbols = 0.6,ylim = c(y_min,y_max), zlim = c(z_min,z_max),asp = -2,
                                      xlab = "sample / dataset",ylab ="X2",zlab = "X1")
  
  for (i in as.numeric(unique(data$sample_id))){
  x0 <- i
  spd$points3d(x=rep(x0,5),y=c(y_min,y_min,y_max,y_max,y_min), z=c(z_min,z_max,z_max,z_min,z_min), type="l", col="black", lwd=0.8)
  }
  
  if(!is.null(pdf_name))dev.off()
}
