## Authors : Pacôme Prompsy
## Title : Wrappers & functions to preprocess & reduce dimensionality of sc Matrix
## Description : Funtions to load, filter, normalize & reduce sc Epigenetic
## Matrices prior to analysis

# Wrapper to create the single cell experiment from sparse datamatrix and annot, remove Zero 
# count Features then Cells and calculate QC Metrics (scater::calculateQCMetrics)
#' Title
#'
#' @param datamatrix 
#' @param annot 
#' @param removeZeroCells 
#' @param removeZeroFeatures 
#'
#' @return
#' @export
#'
#' @examples
create_scExp <- function(datamatrix, annot, removeZeroCells = TRUE, removeZeroFeatures=TRUE){
  
  if(nrow(datamatrix) != nrow(annot))  stop('datamatrix and annot should contain the same number of cells')
  
  scExp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = datamatrix), rowData = annot)
  
  dim_b = dim(scExp)
  if(removeZeroFeatures) scExp <- scExp[,(Matrix::colSums(counts(scExp)>0)>0)] # remove windows that do not have any read in any cells
  if(removeZeroCells) scExp <- scExp[(Matrix::rowSums(counts(scExp)>0)>0), ] # remove cells that do not have any read in any cells
  
  if( dim(scExp)[2] != dim_b[2] ) cat( dim_b[2] - dim(scExp)[2], "features with 0 signals were removed.")
  if( dim(scExp)[1] != dim_b[1] ) cat( dim_b[1] - dim(scExp)[1], "cells with 0 signals were removed.")
  
  scExp <- scater::calculateQCMetrics(scExp)
  
  return(scExp)
}

# Function to filter out cells & features from sparse matrix based on total count per cell,
# number of cells 'ON' (count >= 2) in features and top covered cells that might be doublets
#' Title
#'
#' @param scExp 
#' @param min_cov_cell 
#' @param quant_removal 
#' @param percentMin 
#' @param bin_min_count 
#'
#' @return
#' @export
#'
#' @examples
filter_scExp <- function(scExp, min_cov_cell = 1600, quant_removal = 95, percentMin = 1, bin_min_count = 2){
  
  if(is.null(scExp)) warn("Please specify a SingleCellExperiment")
  
  cellCounts = Matrix::rowSums(counts(scExp))
                               
  thresh <- quantile(cellCounts, probs=seq(0, 1, 0.01))
  
  sel1000 =  ( cellCounts > 1000 & cellCounts < thresh[quant_removal+1] ) 
  sel <- ( cellCounts > min_cov_cell & cellCounts < thresh[quant_removal+1] )
  
  # Create the binary matrix of cells >= 1000 counts
  # if count >= bin_min_count, count = 1 else
  bina_counts <- counts(scExp)[,sel1000]
  sel_below_2  <- (bina_counts < bin_min_count)
  bina_counts[sel_below_2] <- 0
  sel_above_1 = (bina_counts >= bin_min_count)
  bina_counts[sel_above_1] <- 1
  
  nCells_in_feature = Matrix::colSums(bina_counts) 
  fixedFeature <- ( nCells_in_feature > ((percentMin/100)*(nrow(bina_counts))) ) # Feature selection
  
  scExp <- scExp[sel,]
  rowData(scExp) <- rowData(scExp)[sel,]
  
  scExp <- scExp[,fixedFeature]
  
  return(scExp)
}
