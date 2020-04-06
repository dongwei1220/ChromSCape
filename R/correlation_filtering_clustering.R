## Authors : Pacôme Prompsy Title : Wrappers & functions to filter and cluser
## single cell data based on correlation between cells Description : Wrappers &
## functions to filter and cluser single cell data based on correlation between
## cells

correlation_and_hierarchical_clust_scExp <- function(scExp)
  {
  
  stopifnot(is(scExp, "SingleCellExperiment"))
  
  if(is.null(reducedDim(scExp,"PCA"))) stop("Run PCA on the object before correlation") 
    
  pca = reducedDim(scExp,"PCA")
  pca_t <- Matrix::t(pca)
  
  hc_cor = hclust(as.dist(1 - cor(pca_t)), method="ward.D") 
  hc_cor$labels = rep("",length(hc_cor$labels))
  
  scExp@metadata$hc_cor = hc_cor
  reducedDim(scExp,"Cor") <- cor(pca_t)[,hc_cor$order]
  
}

filter_correlated_cell_scExp <- function(scExp, random_iter = 50, corr_threshold = 99, percent_correlation = 1){
  
  stopifnot(is(scExp, "SingleCellExperiment"))
  if(is.null(reducedDim(scExp,"Cor"))) stop("No correlation, run correlation before filtering.")
  
  pca_t = Matrix::t(reducedDim(scExp,"PCA"))
  correlation_values <- vector(length=random_iter)
  corChIP <- reducedDim(scExp,"Cor")
  z <- matrix(sample(pca_t), nrow=dim(pca_t)[1])
  threshold <- quantile(cor(z), probs=seq(0,1,0.01))
  limitC <- threshold[corr_threshold+1]
  
  set.seed(147)
  for(i in 1:random_iter){
    random_mat <-  matrix(sample(pca_t), nrow=dim(pca_t)[1])
    threshold <- quantile(cor(random_mat), probs=seq(0,1,0.01))
    limitC <-  threshold[corr_threshold+1]
    correlation_values[i] = limitC
  }
  
  limitC_mean = mean(correlation_values,na.rm = T)
  
  selection_cor_filtered <- (apply(corChIP, 1, function(x) length(which(x > limitC_mean)))
              ) > (percent_correlation*0.01)*dim(corChIP)[1]
  
  pca_t_cor_filtered <- pca_t[, selection_cor_filtered]
  hc_cor_cor_filtered <- hclust(as.dist(1 - cor(pca_t_cor_filtered)), method="ward.D")
  hc_cor_cor_filtered$labels = rep("",length(hc_cor_cor_filtered$labels))
  cor_cor_filtered <- pca_t_cor_filtered[, hc_cor_cor_filtered$order]
  
  annot_sel <- annot[selection_cor_filtered,]
  
  anocol_sel <- geco.annotToCol2(annotS=annot_sel[, annotCol], annotT=annot_sel, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", dataset_name,"Annotation_legends_reclustering.pdf"), categCol=NULL)
  lapply(colnames(anocol_sel), function(col){ if(col != "total_counts"){ anocol_sel[, col] <<- as.character(reactVal$annotColors[which(reactVal$annotColors$Sample %in% rownames(anocol_sel)), paste0(col, '_Color')]) } })
  
  return(scExp)
}
