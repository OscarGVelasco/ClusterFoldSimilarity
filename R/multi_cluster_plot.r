#' Plot several single-cell datasets/samples dimensionality reductions (PCA, UMAP, tSNE) in a single plot 3D plot, each layer corresponding to a dataset/sample.
#'
#' `multi_cluster_plot()` plot dimensionality reduction (PCA, UMAP, tSNE) of a batch of samples/datasets.
#'
#' This function will plot a list of single-cell datasets/samples dimensionality reduction in a 3D layered plot. 
#'
#' @param sce_list List. A list of single cell experiments. At least 2 single cell experiments are needed.
#' @param similarity_table DataFrame. Optional. A similarity table previously calculated using ClusterFoldSimilarity function.
#' @param dim_method Character. Specify the dimensionality reduction used, either PCA, UMAP or tSNE. Default: PCA.
#' @param pdf_name Character. If specify the plot is saved on a PDF file with the given name. Default: plot in the current graphic device.
#' 
#' @return A dimensionality reduction 3D plot (PCA, UMAP, tSNE) of a batch of samples/datasets.
#' 
#' @import scatterplot3d
#' 
#' @export
multi_cluster_plot <- function(sce_list = NULL,similarity_table = NULL,dim_method="PCA",pdf_name = NULL){

  # sce_list <- pancreas.sc.list
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
    if(dimMeth %in% names(x@reductions)){return(Seurat::Embeddings(x[[dim_method]])[,1:2]) # If the dimensionality has been previously calculated, we obtain the values.
    }else{
        return(useDimFunc(x)) # If the required dimensionality is NOT in the object, we calculate it.
      }}
  }
  if(is_sce){getDimFunc <- function(x,dim_method = dimMeth,useDimFunc=useDimFunc){
    if(dimMeth %in% SingleCellExperiment::reducedDimNames(x)){return(SingleCellExperiment::reducedDim(x,dimMeth)[,1:2]) # If the dimensionality has been previously calculated, we obtain the values.
    }else{
      return(useDimFunc(x)) # If the required dimensionality is NOT in the object, we calculate it.
    }}
  }
  
  # if(is_sce){
  #   if(dimMeth %in% reducedDimNames(sce_list[[1]])){data <- reducedDim(sce_list[[1]],dimMeth)}else{data <- useDimFunc(sce_list[[1]])}
  #   }
  if(is_sce){
    data <- cbind.data.frame(getDimFunc(sce_list[[1]]),sample_id=1, cluster =  SingleCellExperiment::colLabels(sce_list[[1]]))
    for (i in seq_along(sce_list)[-1]){
      data <- rbind.data.frame(data,cbind.data.frame(getDimFunc(sce_list[[i]]),sample_id=i, cluster = SingleCellExperiment::colLabels(sce_list[[i]])))
    }
    }
  if(is_seurat){data <- cbind.data.frame(getDimFunc(sce_list[[1]]),sample_id=1, cluster =  Seurat::colLabels(sce_list[[1]]))}
  
  data$sample_id <- data$sample_id
  y_max <- max(data[,1]) +2
  y_min <- min(data[,1]) -2
  z_max <- max(data[,2]) +2
  z_min <- min(data[,2]) -2
  if(!is.null(pdf_name)){
    graphics.off()
    pdf(file = pdf_name, width = 40, height = 16)
  }
  
  cluster_colors <- c("#D1E5F0","#A6761D","#252525","#0868AC","#08306B","#004529","#0570B0","#2171B5","#A50026","#BCBDDC","#B15928","#8C96C6","#1B7837","#238B45",
                      "#A6BDDB","#DE77AE","#3F007D","#D0D1E6","#313695","#FB6A4A","#B3B3B3","#E08214","#006837","#7BCCC4","#081D58","#F46D43","#78C679","#C7EAE5",
                      "#74C476","#02818A","#67001F","#FDDBC7","#C6DBEF","#E31A1C","#BDBDBD","#FDD49E","#EF6548","#FC8D62","#F7F7F7","#969696","#253494","#D73027",
                      "#FEC44F","#CCEBC5","#053061","#1D91C0","#FB9A99","#FEE0D2","#762A83","#FDBF6F","#D9D9D9","#A8DDB5","#67A9CF","#FEB24C","#B2182B","#FD8D3C",
                      "#800026","#045A8D","#ADDD8E","#542788","#525252")
  
  data <- cbind.data.frame(data, color = "NA")
    
  if(!is.null(similarity_table)){ # If a similarity table is available...  
  
    for (i in seq_along(sce_list)){ # Sample loop - i
    for (j in unique(similarity_table[similarity_table$sample_l==i,]$cluster_l)){ # Cluster loop - j
      # Coloring of the clusters
      if(all(data[(data$sample_id == i) & (data$cluster == j),]$color == "NA")){
        single_color <- sample(cluster_colors,size = 1)
        cluster_colors <- setdiff(cluster_colors,single_color)
        data[data$sample_id == i & data$cluster == j,]$color <- single_color
      }else{
      single_color <- unique(data[(data$sample_id == i) & (data$cluster == j),]$color)
      }
      
      indx <- similarity_table[similarity_table$sample_l==i & similarity_table$cluster_l == j,c("sample_r","cluster_r")]
      
      indx2 <- do.call(rbind.data.frame,
                       apply(indx,1,function(x){
                         similarity_table[(similarity_table$sample_l == x[1]) & (similarity_table$cluster_l == x[2]) &
                                            (similarity_table$sample_r == i) & (similarity_table$cluster_r == j),c("sample_l","cluster_l")]
                       }))
      
      indx3 <- unlist(apply(indx2,1,function(x){
        tmp <- which( (data$sample_id == x[1]) & (data$cluster == x[2]))
        if(all(data[tmp,]$color == "NA")){return(tmp)}
        }))
      
      if(length(indx3) != 0){
        data[indx3,]$color <- single_color
      }
    }}
  } # END of similarity table processing
  
  data_na <- data[data$color == "NA",]
  if(nrow(data_na) != 0){
    # samples_pre <- rownames()
    data_na <- do.call(rbind,unname(by(data_na,data_na[,c("cluster")],function(x){
    single_color <- sample(cluster_colors,size = 1)
    cluster_colors <- setdiff(cluster_colors,single_color)
    x[,"color"] <- single_color
    return(x)
    })))
  data[rownames(data_na),]$color <- data_na$color
  }
  
  spd <- scatterplot3d::scatterplot3d(x = data$sample_id,y = data[,1],z = data[,2],color = data$color,
                                      angle = 60,pch = 20,box = FALSE,col.grid = "grey",
                                      cex.symbols = 0.6,ylim = c(y_min,y_max), zlim = c(z_min,z_max),asp = -2,
                                      xlab = "sample / dataset",ylab ="X2",zlab = "X1")
  
  for (i in as.numeric(unique(data$sample_id))){
  x0 <- i
  spd$points3d(x=rep(x0,5),y=c(y_min,y_min,y_max,y_max,y_min), z=c(z_min,z_max,z_max,z_min,z_min), type="l", col="black", lwd=0.8)
  }
  
  if(!is.null(pdf_name))dev.off()
}
